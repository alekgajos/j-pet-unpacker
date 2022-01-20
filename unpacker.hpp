#ifndef UNPACKER_H
#define UNPACKER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <queue>
#include <fstream>
#include <iostream>

#include "unpacker_types.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-label"

namespace unpacker {

    // flags and constants
    uint32_t verbosity = 0; // 0 - none, 1 - crit. errors, 2 - errors/warnings

    double kCoarsePeriod = 2352.941; // ps
    double kFineBinWidth = kCoarsePeriod / 128; // ps, more or less accurate (active TDC bin count is rather hard to be determined [100-128])
    bool kPerformTDCCalib = true;

    bool kSkipAfterReferenceSample = true; // unpacker filter, rejects samples that arrived after the reference time

    bool kCapTimeAtThreshold = true; // upacker filter, rejects samples that calculated times are larger than kTimeCapThreshold
    uint32_t kTimeCapThreshold = 50000000; // ps

    uint32_t kOverflowLimit = 400; // if there are more samples in packet, overflow is detected and marked as warning in endp_stats_t struct.

    bool kIsOnline = false; // reject additionals headers from online datastream

    // functions
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
            const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib, 
            std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<uint32_t>>> &tdc_calib );

    int32_t get_time_window( 
            meta_t &meta_data,
            std::unordered_map<uint32_t, std::vector<hit_t>> &original_data, 
            std::unordered_map<uint32_t, std::vector<hit_t>> &filtered_data,
            std::unordered_map<uint32_t, std::vector<sigmat_t>> &preproc_data,
            const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib,
            std::ifstream &fp );

    int32_t get_time_window_repaired( 
            meta_t &fixed_meta_data,
            std::unordered_map<uint32_t, std::vector<hit_t>> &fixed_data,
            const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib,
            std::ifstream &fp );

    void calculate_time( 
            uint32_t endp_id,
            std::vector<hit_t> &v, 
            std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<uint32_t>>> &tdc_calib );

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

    static uint32_t byte_cntr = 0;

    unsigned char buffer[4];
    bool ret = true;

    if ( fp.read( (char *)(&buffer[0]), 4) ) {
        *data = buffer[3]<<24 | buffer[2]<<16 | buffer[1]<<8 | buffer[0];
        byte_cntr+=4;
    } else {
        if ( verbosity >  0 ) {
            std::cout << "EOF" << std::endl;
        }

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

    // dump some headers from online analysis
    if ( kIsOnline ) {
        for (uint32_t i=0; i < 7; i++) {
            int32_t succ = unpacker::read_4b( &buff, fp );
            if ( !succ ) {
                return 0;
            }
        }
    }

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
        const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib, 
        std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<uint32_t>>> &tdc_calib ) {

    // loads tdc calibrations from .txt file.
    // input map of <16-bit endpoint id, path to tdc calibration file>, e.g. <0xa110, "/storage/tdc_calibs/0xa110.txt">
    // output 

    for (auto const &pair : paths_to_tdc_calib) {
        uint32_t id = pair.first;
        std::vector<uint32_t> calib;
        std::ifstream calib_file(pair.second);
        
        if (!calib_file) {
            if ( verbosity > 0 ) {
                std::cout << "Error! File " << pair.second.c_str() << " could not be open." << std::endl;
            }

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
    }
    
    return true;
} 


int32_t unpacker::get_time_window( 
        meta_t &meta_data,
        std::unordered_map<uint32_t, std::vector<hit_t>> &original_data, 
        std::unordered_map<uint32_t, std::vector<hit_t>> &filtered_data,
        std::unordered_map<uint32_t, std::vector<sigmat_t>> &preproc_data,
        const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib,
        std::ifstream &fp ) {

    // input: file descriptor pointing to opened hld file,
    // paths to tdc calibrations.
    // output: original, filtered and preproc data if found (by reference), 
    // meta data regarding the time window (by reference),
    // operation status (by value) - 0 in case of EOF, 1 in case of success, -1 in case of failure, 2 in case of error that may be recoverable

    // clear the output structes
    original_data.clear();
    filtered_data.clear();
    preproc_data.clear();
    
    // statics...
    static uint32_t queue_id = 0;
    static bool init_run = true;
    static std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<uint32_t>>> tdc_calib;

    // during first run load the tdc calibration tables
    if ( init_run == true ) {
        init_run = false;

        if ( kPerformTDCCalib ) {
            int32_t succ = load_tdc_calib(paths_to_tdc_calib, tdc_calib);
            
            if ( !succ ) {
                if ( verbosity > 0 ) {
                    std::cout << "Error, TDC mapping tables were not loaded correctly" << std::endl;
                }

                return -1;
            }
        }

        // read file header from online data-stream
        if ( kIsOnline ) {
            uint32_t buff;

            for (uint32_t i=0; i < 16; i++) {
                int32_t succ = unpacker::read_4b( &buff, fp );
                if ( !succ ) {
                    return 0;
                }
            }
        }
    }

    std::vector<uint32_t> raw_data;

    // loop forever till valid queue is found
    loop_till_valid_queue : while (1) {
        int32_t succ = read_queue(raw_data, fp);
        if ( !succ ) { // eof!
            return 0;
        }

        queue_id++;

        if ( raw_data[0] <= 0x20 ) {
            // the queue is too small, contains just the headers of it
            continue;
        } else {
            break;
        }
    }

    // initialize meta data structure
    meta_data.endp_stats.clear();
    meta_data.queue_id = queue_id;
    meta_data.tw_trigger_id = 0;

    // prefill the output maps with raw data
    uint32_t qi = 8; // qi is the queue iterator. It points to the currently
            // preprocessed word. Rest of iterators (such as ci, ei) is
            // just for help. Initially qi is set to 8 to skip first 8 words of queue, 
            // as these are its headers

    // loop over concentrators...
    foreach_concentrator : for (; qi < raw_data.size(); ) {
        // parse 4 words from the concentrator's headers
        for (uint32_t ci=0; ci < 4 && qi < raw_data.size(); ci++, qi++ ) {
            // TODO fill the meta_t. maybe do some comparision
            if ( ci == 2 ) {
                uint32_t buff = reverse_DJ(raw_data[qi]);

                // not a valid concentrator, not a valid endpoint. Probably wandering around...
                if ( (buff & 0xf0ffffff) != 0xA0000000 ) { 
                    break;
                }
            }
            
            if ( ci == 3 ) {
                meta_data.tw_trigger_id = reverse_TID(raw_data[qi]);
                meta_data.raw_trigger_id = reverse_DJ(raw_data[qi]);
            }
        }

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
                continue;
            } else {
                // if here, it means the current data sample does not come from endpoint.
                // it might be concentrator, padding or next queue.
                break; 
            }

            // determine the data type of the packet to unpack it correctly
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
                            
                            // header 1
                            // if ( ei == 0 ) { 
                            // }

                            // header 2
                            if ( ei == 1 ) {
                                // fill the meta data structure
                                meta_data.endp_stats[endp_id].trigger_id = (reverse_TDC(raw_data[qi]) & 0xFF);
                                meta_data.endp_stats[endp_id].overflow_flag = ((reverse_TDC(raw_data[qi]) >> 12) & 0x1);
                                meta_data.endp_stats[endp_id].len = endp_len;
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

                    // run additional stats to see if the reference time is present. Unfortunatelly this needs to be done here.
                    uint32_t ref = get_ref(buff_v);
                    if ( ref == 0 ) {
                        meta_data.endp_stats[endp_id].is_reference = false;
                    } else {
                        meta_data.endp_stats[endp_id].is_reference = true;
                    }

                    calculate_time(endp_id, buff_v, tdc_calib); // fill the time information

                    // load only if the vector size if larger than 0
                    if ( buff_v.size() > 0 ) {
                        original_data[endp_id] = buff_v;
                    }

                    break;

                } // case A 
                
                case 0xB: { // TODO: validate

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

                case 0xC: { // TODO: finish :) and validate

                    std::vector<sigmat_t> buff_v;

                    // there is some sigmat hexes from the endpoint
                    if (endp_len > 2) {

                        uint32_t low_word;
                        uint32_t high_word;

                        foreach_preprocessed_data : for (uint32_t ei=0; qi < raw_data.size() && ei < endp_len; ei++, qi++ ) {
                            // ei is the endpoint word helper-iterator.
                            // qi is the queue iterator. It points to the currently
                            // preprocessed word. Rest of iterators (such as ci, ei) is
                            // just for help. Initially qi is set to 8 to skip first 8 words of queue, 
                            // as these are its headers
                            

                            if (ei != 0 && ei != endp_len - 1) {  // process only data between the header and the trailer

                                if (ei % 2 == 1) { // odd is the lower data part
                                    low_word = reverse_TDC(raw_data[qi]);
                                }
                                else { // even is the upper data part and completes a sigmat
                                    high_word = reverse_TDC(raw_data[qi]);
                                    
                                    sigmat_t buff;
                                    buff.lead_time = low_word & (0x03ffffff);
                                    buff.tot_time = high_word & (0x03ffffff);
                                    buff.strip_id = (low_word & 0x3c000000) >> 26;
                                    buff.multiplicity = (high_word & 0x3c000000) >> 26;
                                    
                                    buff_v.push_back(buff);
                                }
                            }

                            if (buff_v.size() > 0) {
                                preproc_data[endp_id] = buff_v;
                            }
                        }
                    }
                    // no sigmat hexes, skip the headers and proceed
                    else {
                        qi += 2;
                    }

                    break;

                    
                    // foreach_preprocessed_data : for (uint32_t ei=0; qi < raw_data.size() && ei < endp_len; (qi % 2 == 1 ? ei++ : ei=ei), qi++ ) {
                    //     // TODO decoding                                           ^ (!) increment only on even word

                    //     uint32_t low_word;
                    //     uint32_t high_word;

                    //     std::cout<<"w: "<<std::hex<<raw_data[qi]<<std::endl;

                    //     if ( qi % 2 == 0 ) {
                    //         low_word = reverse_TDC(raw_data[qi]);
                    //     } else {
                    //         high_word = reverse_TDC(raw_data[qi]);

                    //         uint64_t tdc_word = ((uint64_t) high_word) << 32 | low_word;

                    //         sigmat_t buff = {
                    //             .lead_time = 0,
                    //             .tot_time = 0,
                    //             .strip_id = 0,
                    //             .multiplicity = 0                               
                    //         };

                    //         std::cout<<"got strip "<<buff.strip_id<<std::endl;

                    //         buff_v.push_back(buff);
                    //     }

                        
                    // } // for

                    // preproc_data[endp_id] = buff_v;

                    // break;

                } // case C

                default: {
                    if ( verbosity > 1 ) {
                        std::cout << "Error, unrecognized sequence" << std::endl;
                    }

                    break;
                    // return -1; // idk, maybe it should return...?
                }

            } // switch

        } // for enpoints

        // skip the padding
        if ( ( qi < raw_data.size() ) && ( raw_data[qi] == 0x05050505 ) ) {
            qi++;
        }

    } // for concentrators

    // run additional check to see if the queue is valid

    // get any trigger_id
    auto i = meta_data.endp_stats.begin();
    uint32_t buff = i->second.trigger_id;

    for ( auto const &pair : meta_data.endp_stats ) {
        if ( buff != pair.second.trigger_id ) {
            return 2;
        }

        buff = pair.second.trigger_id;
    }

    // return '1' if whole queue was read without errors
    return 1;
}


int32_t unpacker::get_time_window_repaired( 
            meta_t &fixed_meta_data,
            std::unordered_map<uint32_t, std::vector<hit_t>> &fixed_data,
            const std::unordered_map<uint32_t, std::string> &paths_to_tdc_calib,
            std::ifstream &fp ) {

    // input: file descriptor pointing to opened hld file,
    // paths to tdc calibrations.
    // output: fixed data if found (by reference), 
    // meta data regarding the time window (by reference),
    // operation status (by value) - 0 in case of EOF, 1 in case of success, -1 in case of failure, 2 if the fix was performed successfully.

    static std::unordered_map<uint32_t, std::vector<hit_t>> prev_original_data;
    static meta_t prev_meta_data;

    std::unordered_map<uint32_t, std::vector<hit_t>> original_data;
    meta_t meta_data;

    // dummy variables required by 'get_time_window' function
    std::unordered_map<uint32_t, std::vector<hit_t>> filtered_data;
    std::unordered_map<uint32_t, std::vector<sigmat_t>> preproc_data;

    // clear the input structs
    fixed_meta_data.endp_stats.clear();
    fixed_data.clear();

    int32_t ret = get_time_window(meta_data, original_data, filtered_data, preproc_data, paths_to_tdc_calib, fp);

    if ( ret == -1 || ret == 0 ) { // return normally on error or eof. There's nothing to be fixed.
        return ret;
    }

    if ( ret == 1 ) {
        prev_original_data.clear();
        prev_original_data = original_data;
        prev_meta_data.endp_stats.clear();
        prev_meta_data = meta_data;

        // assing outputs
        fixed_data = original_data;
        fixed_meta_data = meta_data;

        return ret;
    }
    
    // at this point the 'ret' from 'get_time_window' was 2, try to fix

    // flag of the quality of the fix.
    bool is_fixed = true;

    if ( ret == 2 ) { 
        std::unordered_map<uint32_t, uint32_t> late_endpoints; // those endpoints arrived late. The tw of the rest of endpoints must be taken from previous tw
        std::unordered_map<uint32_t, uint32_t> fixed_endpoints; // will hold the list of endpoints after the fix. Used to check the fix quality
                // ^ endp_id, ^ trigger_id
        
        // for endpoints that were late...
        for ( auto const &pair : meta_data.endp_stats ) {
            if ( pair.second.trigger_id != ((meta_data.tw_trigger_id) & 0xFF) ) { // endpoint tw is a mismatch (endpoints was late)
                late_endpoints[pair.first] = pair.second.trigger_id; // append a list of endpoints that are from current tw.

                // save the data of this endpoint
                fixed_endpoints[pair.first] = pair.second.trigger_id;
                fixed_meta_data.endp_stats[pair.first] = pair.second;

                if ( original_data.count(pair.first) ) { // make sure there's anything to be assigned from
                    fixed_data[pair.first] = original_data[pair.first];
                }
            }
        }

        // for the rest of the endpoints
        for ( auto const &pair : prev_meta_data.endp_stats ) { // take data from previous tw
            if ( !late_endpoints.count(pair.first) ) { // must not be in 'late_endpoints' list
                // save the data of this endpoint
                fixed_endpoints[pair.first] = pair.second.trigger_id;
                fixed_meta_data.endp_stats[pair.first] = pair.second;

                if ( prev_original_data.count(pair.first) ) { // make sure there's anything to be assigned from
                    fixed_data[pair.first] = prev_original_data[pair.first];
                }
            }
        }

        // validate the fix...       
        if ( fixed_endpoints.size() != 0 ) {
            // get any trigger_id in 'fixed_endpoints' list
            auto i = fixed_endpoints.begin();
            uint32_t buff = i->second; 

            for ( auto const &pair : fixed_endpoints ) {
                if ( buff != pair.second ) { // at least one trigger_id in 'fixed_endpoints' list differ...
                    if ( verbosity > 1 ) {
                        std::cout << "Error, time-window " << std::hex << meta_data.tw_trigger_id << " could not be fixed..." << std::endl;
                    }
                    // mark fixed as invalid
                    is_fixed = false;
                }
            }
        } else {
            is_fixed = false;
        }

        // assign the rest of meta_data infos
        fixed_meta_data.queue_id = meta_data.queue_id;
        fixed_meta_data.tw_trigger_id = prev_meta_data.tw_trigger_id; // ... valid?
    }

    // fill the 'previous window' data structes
    prev_original_data.clear();
    prev_original_data = original_data;
    prev_meta_data.endp_stats.clear();
    prev_meta_data = meta_data;

    // at this point the ret from the 'get_time_window' can be only '2'.
    // return '2' in case of successful fix, '-1' on fail (time-window should not be used)

    return is_fixed ? 2 : -1;
   
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
        std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<uint32_t>>> &tdc_calib ) {

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

#pragma GCC diagnostic pop
#endif
