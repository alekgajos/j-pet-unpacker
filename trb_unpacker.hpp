#ifndef UNPACKER_H
#define UNPACKER_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <vector>

#include "tdc_calib.hpp"
#include "unpacker_types.hpp"

namespace trb_unpacker {

size_t kWordSize = 4;

typedef struct EventHdr {
  uint32_t fullSize;
  uint32_t decoding;
  uint32_t id;
  uint32_t seqNr;
  uint32_t date;
  uint32_t time;
  uint32_t runNr;
  uint32_t pad;
} hdr_t;

typedef struct SubEventHdr {
  uint32_t size;
  uint32_t decoding;
  uint32_t hubAddress;
  uint32_t trgNr;
} subhdr_t;

hdr_t event_header;
subhdr_t subevent_header;

// functions
inline size_t align8(size_t length) { return 8 * size_t((length - 1) / 8 + 1); }

uint32_t reverse_word(uint32_t word);

size_t process_single_tdc(
    uint32_t *data, size_t data_size,
    std::unordered_map<uint32_t, std::vector<unpacker::hit_t>> &original_data);

unpacker::hit_t
calculate_time(uint32_t word, unsigned int channel_number, int epoch,
               double reference_time,
               std::unordered_map<uint32_t, std::vector<uint32_t>> &tdc_calib);

int32_t get_time_window(
    unpacker::meta_t &meta_data,
    std::unordered_map<uint32_t, std::vector<unpacker::hit_t>> &original_data,
    std::unordered_map<uint32_t, std::vector<unpacker::hit_t>> &filtered_data,
    std::unordered_map<uint32_t, std::vector<unpacker::sigmat_t>> &preproc_data,
    std::istream &fp);

} // namespace trb_unpacker

// returns: nunber of words used from the data
size_t trb_unpacker::process_single_tdc(
    uint32_t *data, size_t data_size,
    std::unordered_map<uint32_t, std::vector<unpacker::hit_t>> &original_data) {

  size_t words_used = 0;
  uint32_t *word_it = data;
  std::vector<unpacker::hit_t> hits;

  uint32_t word = reverse_word((*word_it));
  words_used++;
  word_it++; // TODO: use post-increment above

  uint32_t tdc_number = word & 0xffff;
  std::cout << "TDC number = " << tdc_number << std::endl;
  int32_t internal_size = word >> 16;
  std::cout << "Internal size = " << internal_size << std::endl;

  if (unpacker::tdc_calib::perform_tdc_calib &&
      !unpacker::tdc_calib::tdc_calib.count(tdc_number)) {
    // skip ahead
    words_used += internal_size + 1;
    word_it += internal_size + 1;
    return words_used;
  }

  int epoch = -1;
  int channel, prev_channel = -1;
  uint32_t coarse, prev_coarse;
  uint32_t fine, prev_fine;
  int repetition_counter = 0;
  double reference_time = 0.;

  while (internal_size-- >= 0) {

    word = reverse_word((*word_it));
    words_used++;
    word_it++; // TODO: user post-increment above

    uint32_t header = word >> 29;

    switch (header) {
    case 3: // epoch counter
      epoch = word & 0xffffffff;
      break;
    case 4: // time data
    {
      channel = (word >> 22) & 0x7f;

      coarse = word & 0x7ff;
      fine = (word >> 12) & 0x3ff;

      // count possible repetitions as a sign of corrupted data
      if (channel > 0) {
        if (fine == prev_fine && coarse == prev_coarse &&
            channel == prev_channel) {
          repetition_counter++;
        }
        // TODO: handle max number of allowed repetitions
      }

      prev_channel = channel;
      prev_coarse = coarse;
      prev_fine = fine;

      unpacker::hit_t hit =
          calculate_time(word, channel, epoch, reference_time,
                         unpacker::tdc_calib::tdc_calib[tdc_number]);
      // handle reference
      if (channel == 0) {
        reference_time = hit.time;
      }

      hits.push_back(hit);
      break;
    }
    default:
      break;
    }
  }

  // keep only if the vector size if larger than 0
  if (hits.size() > 0) {
    original_data[tdc_number] = hits;
  }

  return words_used;
}

int32_t trb_unpacker::get_time_window(
    unpacker::meta_t &meta_data,
    std::unordered_map<uint32_t, std::vector<unpacker::hit_t>> &original_data,
    std::unordered_map<uint32_t, std::vector<unpacker::hit_t>> &filtered_data,
    std::unordered_map<uint32_t, std::vector<unpacker::sigmat_t>> &preproc_data,
    std::istream &fp) {

  // input: file descriptor pointing to opened hld file,
  // output: original data (by reference), filtered and preproc data is always
  // empty in case of trb daq meta data regarding the time window (by reference)
  // - always empty in case of trb daq, operation status (by value) - 0 in case
  // of EOF, 1 in case of success, -1 in case of failure

  // clear the output structs
  original_data.clear();
  filtered_data.clear();
  preproc_data.clear();

  // statics...
  static uint32_t queue_id = 0;
  static bool init_run = true;

  // during first run load the tdc calibration tables
  if (init_run == true) {
    init_run = false;

    // skip the file header
    fp.ignore(32);
  }

  std::vector<uint32_t> raw_data;
  int event_size = 0;

  // loop forever till valid non-empty event is found
  while (true) {

    if (!fp.read((char *)&event_header, sizeof(hdr_t))) {
      return 0;
    }

    queue_id++;
    event_size = event_header.fullSize;

    if (event_size == 32) { // the event is header-only
      continue;
    } else {
      break;
    }
  }
  std::cout << "Event size: " << event_size << std::endl;

  int remaining_bytes = event_size - sizeof(hdr_t);

  while (remaining_bytes > 0) { // loop over data modules

    if (!fp.read((char *)&subevent_header, sizeof(subhdr_t))) {
      return 0;
    }

    int data_size = reverse_word(subevent_header.size) - sizeof(subhdr_t);
    uint32_t padded_data_size = align8(data_size);
    uint32_t padding_bytes = padded_data_size - data_size;
    std::cout << "Data size: " << data_size << std::endl;
    std::cout << "Padded size: " << padded_data_size << std::endl;
    uint32_t *data =
        new uint32_t[data_size]; // TODO: reuse once-allocated memory

    if (!fp.read((char *)(data), data_size)) {
      return 0;
    }
    uint32_t *data_it = data;

    if (*data_it != 0) {

      uint32_t module_id = reverse_word(subevent_header.hubAddress);
      std::cout << "Address: " << module_id << std::endl;
      std::vector<unpacker::hit_t> hits;
      size_t processed_words = 0;

      // process TDC data
      while (data_size > 0) {
        processed_words = process_single_tdc(data_it, data_size, original_data);
        data_size -= kWordSize * processed_words;
        data_it += processed_words;
      }

    } else {
      std::cerr << "First data word is empty, skipping event." << std::endl;
    }

    // done processing TDC data
    delete[] data;
    fp.ignore(padding_bytes);
    remaining_bytes -= padded_data_size;
    remaining_bytes -= sizeof(subhdr_t);

    std::cout << "Remains: " << remaining_bytes << std::endl;
  }

  return 1; // success
}

unpacker::hit_t trb_unpacker::calculate_time(
    uint32_t word, unsigned int channel_number, int epoch,
    double reference_time,
    std::unordered_map<uint32_t, std::vector<uint32_t>> &tdc_calib) {

  // calculate times
  bool is_rising = (word >> 11) & 0x1;
  int coarse = word & 0x7ff;
  unsigned int fine = (word >> 12) & 0x3ff;
  double full_time = 0;

  std::cout << "Fine time = " << fine << std::endl;

  // TODO: fetch TDC correction
  if (unpacker::tdc_calib::perform_tdc_calib &&
      tdc_calib[channel_number].size() > 0) {
    assert(tdc_calib[channel_number].size() > fine + 1);
    fine = tdc_calib[channel_number][fine + 1];
  } else {
    fine *= 10;
  }

  if (fine == 0x3ff) {
    return unpacker::hit_t{0, 0., -1, -1};
  }

  full_time = (double)(((epoch << 11) * 5.0));
  full_time += ((coarse * 5000.) - fine) / 1000.;
  full_time -= reference_time;

  unpacker::hit_t decoded_hit;
  decoded_hit.time = full_time;
  decoded_hit.is_falling_edge = !is_rising;
  decoded_hit.channel_id = channel_number;

  return decoded_hit;
}

inline uint32_t trb_unpacker::reverse_word(uint32_t word) {
  uint32_t a, b, c, d, e;
  a = word & 0x000000ff;
  b = word & 0x0000ff00;
  c = word & 0x00ff0000;
  d = word & 0xff000000;

  a <<= 24;
  b <<= 8;
  c >>= 8;
  d >>= 24;

  e = a | b | c | d;

  return e;
}

#endif
