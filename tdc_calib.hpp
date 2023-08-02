#ifndef TDC_CALIB_H
#define TDC_CALIB_H

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace unpacker {
namespace tdc_calib {

typedef std::unordered_map<uint32_t,
                           std::unordered_map<uint32_t, std::vector<uint32_t>>>
    tdc_calib_t;

tdc_calib_t tdc_calib;
bool perform_tdc_calib = false;

bool load_tdc_calib(
    const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib);
inline void set_tdc_calib(const tdc_calib_t &calib_data);

} // namespace tdc_calib
} // namespace unpacker

inline void unpacker::tdc_calib::set_tdc_calib(const tdc_calib_t &calib_data) {
  tdc_calib = calib_data;
  perform_tdc_calib = true;
}

bool unpacker::tdc_calib::load_tdc_calib(
    const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib) {

  // loads tdc calibrations from .txt files.
  // input map of <16-bit endpoint id, path to tdc calibration file>, e.g.
  // <0xa110, "/storage/tdc_calibs/0xa110.txt">

  for (auto const &pair : paths_to_tdc_calib) {
    uint32_t id = pair.first;
    std::vector<uint32_t> calib;
    std::ifstream calib_file(pair.second);

    if (!calib_file) {
        std::cout << "Error! File " << pair.second.c_str()
                  << " could not be open." << std::endl;

      perform_tdc_calib = false;
      return false;
    }

    uint32_t buff;
    uint32_t loop_cntr = 0;
    while (calib_file >> buff) {
      loop_cntr++;
      calib.push_back(buff);

      if (loop_cntr % 128 == 0) {
        tdc_calib[id][loop_cntr / 128 - 1] = calib;
        calib.clear();
      }
    }
  }

  perform_tdc_calib = true;
  return true;
}

#endif
