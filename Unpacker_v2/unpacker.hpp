#ifndef UNPACKER_H
#define UNPACKER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

namespace unpacker {

    double kCoarsePeriod = 2352.941; // ps
    double kFineBinWidth = kCoarsePeriod / 128; // ps, more or less accurate
    bool kPerformTDCCalib = true;

    bool kSkipAfterReferenceSample = true; // unpacker filter, rejects samples that arrived after the reference time

    bool kCapTimeAtThreshold = true; // upacker filter, rejects samples that calculated time are larger than kTimeCapThreshold
    uint32_t kTimeCapThreshold = 50000000; // ps

    typedef struct {
        // TODO: put more debug information
        uint32_t trigger_id;
    } meta_t;

    typedef struct {
        uint32_t sample; /* original sample from hld file */
        double time;
        int32_t channel_id;
        int32_t is_falling_edge;
    } hit_t;

    typedef struct {
        double lead_time;
        double tot_time;
        int32_t strip_id;
        int32_t multiplicity;
    } sigmat_t;

    inline uint32_t reverse_TDC( uint32_t sample );
    inline uint32_t reverse_TID( uint32_t sample );
    inline uint32_t reverse_DJ( uint32_t sample );
    inline int32_t get_channel( uint32_t sample );
    inline uint32_t get_coarse( uint32_t sample );
    inline uint32_t get_fine( uint32_t sample );
    inline int32_t get_edge( uint32_t sample );

    inline bool read_4b( uint32_t *data, std::ifstream &fp );
    int32_t read_queue( std::vector<uint32_t> &data, std::ifstream &fp );

    bool load_tdc_calib(
            const std::map<uint32_t, std::string> &paths_to_tdc_calib, 
            std::map<uint32_t, std::map<uint32_t, std::vector<uint32_t>>> &tdc_calib );

    inline int32_t get_time_window( 
            meta_t &meta_data,
            std::map<uint32_t, std::vector<hit_t>> &original_data, 
            std::map<uint32_t, std::vector<hit_t>> &filtered_data,
            std::map<uint32_t, std::vector<sigmat_t>> &preproc_data,
            const std::map<uint32_t, std::string> &paths_to_tdc_calib,
            std::ifstream &fp );

    void calculate_time( 
            uint32_t endp_id,
            std::vector<hit_t> &v, 
            std::map<uint32_t, std::map<uint32_t, std::vector<uint32_t>>> &tdc_calib );

    uint32_t get_ref( const std::vector<hit_t> &v );

}; // namespace unpacker


inline uint32_t unpacker::reverse_DJ( uint32_t sample ) {
    // input: 4-byte hld-word from DJ region (e.g. headers of queue) of hld file;
    // output: 4-byte decoded word
    
    uint32_t a, b, c, d, e;

    a = sample & 0xFF;
    b = sample & 0xFF00;
    c = sample & 0xFF0000;
    d = sample & 0xFF000000;

    a <<= 8;
    b >>= 8;
    c <<= 8;
    d >>= 8;

    e = a|b|c|d;

    return e;
}


inline uint32_t unpacker::reverse_TID( uint32_t sample ) {
    // input: 4-byte hld-word encoding the trigger ID;
    // output: 4-byte decoded word

    uint32_t a, b, c, d, e;

    a = sample & 0xFF;
    b = sample & 0xFF00;
    c = sample & 0xFF0000;
    d = sample & 0xFF000000;

    a <<= 24;
    b <<= 8;
    c >>= 8;
    d >>= 16;

    e = (a|b|c) >> 8;

    return e;
}


inline uint32_t unpacker::reverse_TDC( uint32_t sample ) {
    // input: 4-byte hld-word from TDC region (e.g. data, headers of ftabs) of hld file;
    // output: 4-byte decoded word

    uint32_t a, b, c, d, e;

    a = sample & 0xFF;
    b = sample & 0xFF00;
    c = sample & 0xFF0000;
    d = sample & 0xFF000000;

    a <<= 24;
    b <<= 8;
    c >>= 8;
    d >>= 24;

    e = a|b|c|d;
    
    return e;
}

inline int32_t unpacker::get_channel( uint32_t sample ) { return ((sample >> 24) & 0x7F); }

inline uint32_t unpacker::get_coarse( uint32_t sample ) { return ((sample >> 7) & 0xFFFF); }

inline uint32_t unpacker::get_fine( uint32_t sample ) { return (sample & 0x7F); }

inline int32_t unpacker::get_edge( uint32_t sample ) { return (sample >> 31); }

inline bool unpacker::read_4b( uint32_t *data, std::ifstream &fp ) {
    // input: file descriptor pointing to opened hld file.
    // output: operation status (by value) - 0 in case of EOF, 1 in case of success;
    // 4-byte hld word (by pointer)

    unsigned char buffer[4];
    bool ret = true;

    if ( fp.read( (char *)(&buffer[0]), 4) ) {
        *data = buffer[3]<<24 | buffer[2]<<16 | buffer[1]<<8 | buffer[0];
    } else {
        std::cout << "EOF" << std::endl;
        ret = false;
    }

    return ret;
}


int32_t unpacker::read_queue( std::vector<uint32_t> &data, std::ifstream &fp ) {
    // input: file descriptor pointing to opened hld file
    // output: std::vector filled with raw data from single EventBuilder queue (by reference),
    // operation status (by value) - 0 eof, 1 success

    uint32_t buff = 0;

    // make sure vector is empty at the begining
    data.clear();

    // read first queue word which contains queue's length
    int32_t succ = unpacker::read_4b( &buff, fp );
    if ( !succ ) {
        return 0;
    }

    data.push_back(buff);
    int32_t remaining_bytes = buff - 4;

    // read the rest of the queue
    while ( remaining_bytes > 0 ) {
        succ = unpacker::read_4b( &buff, fp );

        if ( !succ ) {
            return 0;
        }

        data.push_back(buff);
        remaining_bytes -= 4;
    }

    return 1;
}


bool unpacker::load_tdc_calib( 
        const std::map<uint32_t, std::string> &paths_to_tdc_calib, 
        std::map<uint32_t, std::map<uint32_t, std::vector<uint32_t>>> &tdc_calib ) {

    // loads tdc calibrations from .txt file.
    // input map of <16-bit endpoint id, path to tdc calibration file>, e.g. <0xa110, "/storage/tdc_calibs/0xa110.txt">
    // output 

    for (auto const &pair : paths_to_tdc_calib) {
        uint32_t id = pair.first;
        std::vector<uint32_t> calib;
        std::ifstream calib_file(pair.second);
        
        if (!calib_file) {
            std::cout << "Error! File " << pair.second.c_str() << " could not be open." << std::endl;
            return 0;
        }

        uint32_t buff;
        uint32_t loop_cntr = 0;
        while (calib_file >> buff) {
            loop_cntr++;
            calib.push_back(buff);

            if ( loop_cntr % 128 == 0 ) {
                tdc_calib[id][loop_cntr/128 - 1] = calib;
                calib.clear();
            } 
        }

        // if ( loop_cntr != 128*105 ) { // return on error
        //     return 0;
        // }
    }
    
    return true;
} 


int32_t unpacker::get_time_window( 
        meta_t &meta_data,
        std::map<uint32_t, std::vector<hit_t>> &original_data, 
        std::map<uint32_t, std::vector<hit_t>> &filtered_data,
        std::map<uint32_t, std::vector<sigmat_t>> &preproc_data,
        const std::map<uint32_t, std::string> &paths_to_tdc_calib,
        std::ifstream &fp ) {

    // input: file descriptor pointing to opened hld file,
    // paths to tdc calibrations.
    // output: original, filtered and preproc data if found (by reference), 
    // meta data regarding the time window (by reference),
    // operation status (by value) - 0 in case of EOF, 1 in case of success, -1 in case of failure.

    static bool init_run = true;
    static std::map<uint32_t, std::map<uint32_t, std::vector<uint32_t>>> tdc_calib;

    // during first run load the tdc calibration tables
    if ( init_run == true ) {
        init_run = false;
        int32_t succ = load_tdc_calib(paths_to_tdc_calib, tdc_calib);
        
        if ( !succ ) {
            std::cout << "Error, TDC mapping tables were not loaded correctly" << std::endl;
            return -1;
        }
    }

    std::vector<uint32_t> raw_data;

    // loop forever till valid queue is found
    loop_till_valid_queue : while (1) {
        int32_t succ = read_queue(raw_data, fp);
        if ( !succ ) { // eof!
            return 0;
        }

        if ( raw_data[0] <= 0x20 ) {
            // the queue is too small, contains just the headers of it
            continue;
        } else {
            break;
        }
    }

    // int i = 0;
    // for (auto const &val : raw_data) {
    //     std::cout << std::dec << i << " " << std::hex << reverse_TDC(val) << std::endl;
    //     i++;
    // }

    // prefill the output maps with raw data
    uint32_t qi = 8; // qi is the queue iterator. It points to the currently
            // preprocessed word. Rest of iterators (such as ci, ei) is
            // just for help. Initially qi is set to 8 to skip first 8 words of queue, 
            // as these are its headers

    uint32_t error_cntr = 0;

    // loop over concentrators...
    foreach_concentrator : for (; qi < raw_data.size(); ) {
        // parse 4 words from the concentrator's headers
        for (uint32_t ci=0; ci < 4 && qi < raw_data.size(); ci++, qi++ ) {
            // TODO fill the meta_t. maybe do some comparision
            if ( ci == 2 ) {
                uint32_t buff = reverse_DJ(raw_data[qi]);

                // not a valid concentrator, not a valid endpoint. Probably wandering around...
                if ( (buff & 0xf0ffffff) != 0xA0030000 ) { 
                    break;
                }
            }

            if ( ci == 3 ) {
                meta_data.trigger_id = reverse_TID(raw_data[qi]);
            }
        }

        uint32_t endp_trg_num = 0;

        // loop over endpoints...
        foreach_endpoint : for (; qi < raw_data.size(); ) {
            // get endpoint data size (read the first word of the endpoint data).
            // NOTE: after this read 'qi' is not incremented, so this word is being read twice:
            // here and in 'foreach_filtered_data/foreach_original_data/foreach_preprocessed_data' loop.

            uint32_t buff = reverse_TDC(raw_data[qi]);
            uint32_t endp_id = buff >> 16;
            uint32_t endp_len = buff & 0xffff;

            // before moving on check if this really is endpoint
            if ( (endp_id & 0xF00F) == 0xA000 || (endp_id & 0xF00F) == 0xB000 || (endp_id & 0xF00F) == 0xC000) {
                // do nothing, this is endpoint
            } else if ( endp_id == 0xffff && endp_len != 0xffff ) {
                qi += endp_len; // skip the invalid name endpoint
            } else {
                // if here, it means the current data sample does not come from endpoint.
                // it might be concentrator, padding or next queue.
                break; 
            }

            // determin the data type of the packet to unpack it correctly
            uint32_t data_type = endp_id >> 12; // will result A (original), B (prefiltered), C (preprocessed)

            switch (data_type) {

                case 0xA: {

                    std::vector<hit_t> buff_v;

                    // loop over data in endpoint...
                    foreach_original_data : for (uint32_t ei=0; qi < raw_data.size() && ei < endp_len; ei++, qi++ ) {
                        // ei is the endpoint word helper-iterator.
                    
                        if ( ei < 2 || ei >= endp_len - 2 ) {
                            // skip the two words of headers, and two
                            // words of trailers in the endpoint.
                            
                            //
                            if ( ei == 1 ) {
                                if ( endp_trg_num != 0 && endp_trg_num != reverse_TDC(raw_data[qi]) ) {
                                    std::cout << "Error, trigger mismatch on endpoint " << std::hex << endp_id << 
                                            ". Trigger should be " << endp_trg_num << ", but is " << reverse_TDC(raw_data[qi]) << std::endl;
                                    return -1;
                                }
                                endp_trg_num = reverse_TDC(raw_data[qi]);
                            }

                            continue; 
                        }

                        uint32_t tdc_word = reverse_TDC(raw_data[qi]);

                        hit_t buff = {
                            .sample = tdc_word,
                            .time = -1, // initialize it to '-1', as it cannot be calculated now.
                            .channel_id = get_channel(tdc_word),
                            .is_falling_edge = get_edge(tdc_word)  
                        };
                        
                        buff_v.push_back(buff);
                    } // for

                    calculate_time(endp_id, buff_v, tdc_calib); // fill the time information

                    // load only if the vector size if larger than 0
                    if ( buff_v.size() > 0 ) {
                        original_data[endp_id] = buff_v;
                    }

                    break;

                } // case A
                
                case 0xB: {

                    std::vector<hit_t> buff_v;

                    foreach_filtered_data : for (uint32_t ei=0; qi < raw_data.size() && ei < endp_len; ei++, qi++ ) {

                        uint32_t tdc_word = reverse_TDC(raw_data[qi]);

                        hit_t buff = {
                            .sample = tdc_word,
                            .time = -1,
                            .channel_id = get_channel(tdc_word),
                            .is_falling_edge = get_edge(tdc_word)
                        };
                        
                        buff_v.push_back(buff);
                    } // for 

                    calculate_time(endp_id, buff_v, tdc_calib); // fill the time information
                    filtered_data[endp_id] = buff_v;

                    break;

                } // case B

                case 0xC: {

                    std::vector<sigmat_t> buff_v;

                    foreach_preprocessed_data : for (uint32_t ei=0; qi < raw_data.size() && ei < endp_len; (qi % 2 == 1 ? ei++ : ei=ei), qi++ ) {
                        // TODO decoding                                           ^ (!) increment only on even word

                        uint32_t low_word;
                        uint32_t high_word;

                        if ( qi % 2 == 0 ) {
                            low_word = reverse_TDC(raw_data[qi]);
                        } else {
                            high_word = reverse_TDC(raw_data[qi]);

                            uint64_t tdc_word = ((uint64_t) high_word) << 32 | low_word;

                            sigmat_t buff = {
                                .lead_time = 0,
                                .tot_time = 0,
                                .strip_id = 0,
                                .multiplicity = 0                               
                            };

                            buff_v.push_back(buff);
                        }

                        
                    } // for

                    preproc_data[endp_id] = buff_v;

                    break;

                } // case C

                default: {
                    // std::cout << "Warning: endpoint id " << std::hex << endp_id << " of length " << endp_len << " not recognized." << std::endl;
                    break;
                    // return -1; // idk, maybe it should return...?
                }

            } // switch

        } // for enpoints

        // skip the padding
        if ( raw_data[qi] == 0x05050505 ) {
            qi++;
        }

    } // for concentrators

    // return '1' if whole queue was read without errors
    return 1;
}


uint32_t unpacker::get_ref( const std::vector<hit_t> &v ) {

    // input: vector of tdc samples,
    // output: reference channel
    
    if ( v.size() > 0 && get_channel( v[v.size()-1].sample ) == 104 ) {
        return v[v.size()-1].sample; // usually, the reference sample should be on the last position in vector
    }

    // but if it's not true - iterate through the entire vector
    for (auto const &elem : v) {
        if ( get_channel(elem.sample) == 104 ) {
            return elem.sample;
        }
    }

    return 0;
}


void unpacker::calculate_time( 
        uint32_t endp_id,
        std::vector<hit_t> &v, 
        std::map<uint32_t, std::map<uint32_t, std::vector<uint32_t>>> &tdc_calib ) {

    // input: endpoint id, vector of tdc samples, tdc-calibration maps.
    // output: modifies the input vector. Fills it with the information
    // of time.

    // delete the endpoint data if its ID is not present in tdc calibrations.
    if ( !tdc_calib.count(endp_id) ) {
        v.clear();
        return;
    }

    uint32_t ref = get_ref(v);
    
    // Reference sample's edge must be "leading".
    // Calculating times from falling edge sample will result
    // in incorrect values (mostly larger than 50 [us]).

    if ( ref == 0 || get_edge(ref) == 1 ) {
        v.clear();
        return;
    }

    double r, t;

    // calculate the counter values for the reference sample
    if ( kPerformTDCCalib ) {
        r = get_coarse(ref) * kCoarsePeriod - tdc_calib[endp_id][get_channel(ref)][get_fine(ref)];
    } else {
        r = get_coarse(ref) * kCoarsePeriod - get_fine(ref) * kFineBinWidth;
    }

    auto i = std::begin(v);
    bool invalid_samples = false;

    while (i != std::end(v)) {  
        // erase in case of coding error
        if ( i->channel_id > 104 || (get_fine(i->sample) == 127 && get_channel(i->sample) != 104) ) {

            i = v.erase(i);

        } else {
            
            // erase sample in case of incorrect occurrence in a packet, i.e. sample came after reference time
            if ( kSkipAfterReferenceSample ) {
                if ( invalid_samples ) {
                    i = v.erase(i);
                    continue;
                }
                
                if ( i->channel_id == 104 ) {
                    invalid_samples = true;
                }
            }

            // calculate the counter value for the sample depending on calibrations used or not
            if ( kPerformTDCCalib ) {
                t = get_coarse(i->sample) * kCoarsePeriod - tdc_calib[endp_id][(i->channel_id)][get_fine(i->sample)];
            } else {
                t = get_coarse(i->sample) * kCoarsePeriod - get_fine(i->sample) * kFineBinWidth;
            }

            // calculate the time of the hit, take care of counters overflow
            if ( r < t ) {
                i->time = r + (0x10000 * kCoarsePeriod) - t;
            } else {
                i->time = r-t;
            }

            if ( kCapTimeAtThreshold && i->channel_id != 104 ) { // don't delete reference sample
                if ( i->time > kTimeCapThreshold ) {
                    i = v.erase(i);
                    continue;
                }
                
            }

            // for the reference channel overwrite the time with counters value
            if ( i->channel_id == 104 ) {
                i->time = r;
            }

            i++;

        }
    }
}

#endif