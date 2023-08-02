#include "trb_unpacker.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <cstdio>

int main(int argc, char **argv)
{
    std::ifstream fp( argv[1], std::ios::in | std::ios::binary );

    unpacker::meta_t meta_data;
    std::unordered_map<unsigned int, std::vector<unpacker::hit_t>> original_data;
    std::unordered_map<unsigned int, std::vector<unpacker::hit_t>> filtered_data;
    std::unordered_map<unsigned int, std::vector<unpacker::sigmat_t>> preproc_data;
    std::unordered_map<unsigned int, std::string> paths_to_tdc_calib;

    int succ = 1;
    while( succ ) {
        succ = trb_unpacker::get_time_window( meta_data, original_data, filtered_data, preproc_data, fp );

        /* print original data */
         for (auto const &pair: original_data)
         {
             printf("{%02x\n", pair.first);

             for (auto const &val : pair.second)
             {
                 printf("\t%02x (%d, %d, %.0f),\n", val.sample,
                                                   val.is_falling_edge,
                                                   val.channel_id,
                                                   val.time);
             }
          
             printf("}\n");
         }
    
    }


    return 0;
}
