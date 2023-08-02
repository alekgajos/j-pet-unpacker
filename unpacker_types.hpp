#ifndef UNPACKER_TYPES_H
#define UNPACKER_TYPES_H

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace unpacker {


    typedef struct {
        uint32_t trigger_id;
        bool overflow_flag;
        bool is_reference;
        uint32_t len;
    } endp_stats_t;

    typedef struct {
        // basic info
        uint32_t queue_id;
        uint32_t tw_trigger_id;

        // advanced
        std::unordered_map<uint32_t, endp_stats_t> endp_stats;
 
        // temp
        uint32_t raw_trigger_id;
    } meta_t;

    typedef struct {
        uint32_t sample; // original word from hld file
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
  
} // namespace unpacker

#endif
