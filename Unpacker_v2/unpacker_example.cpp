#include "unpacker.hpp"
#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char **argv)
{
    std::ifstream fp( argv[1], std::ios::in | std::ios::binary );

    std::vector<unsigned int> data;

    unpacker::meta_t meta_data;
    std::map<unsigned int, std::vector<unpacker::hit_t>> original_data;
    std::map<unsigned int, std::vector<unpacker::hit_t>> filtered_data;
    std::map<unsigned int, std::vector<unpacker::sigmat_t>> preproc_data;
    std::map<unsigned int, std::string> paths_to_tdc_calib;
    
    for (int i=1; i<=4; i++)
        for (int j=1; j<=12; j++)
        {
            unsigned int id = i << 4 | j;
            std::stringstream ss;
            ss << std::hex << id;
            std::string s_id;
            ss >> s_id;    
            // std::string path_to_calib = "/storage/tdc_calib/A" + s_id + "0_tdc_calib.txt";
            std::string path_to_calib = "/home/baciek/Documents/djpet_software/pyhld/calib/0xa" + s_id + "0_calib.txt";
            paths_to_tdc_calib[ 0xa<<12 | i<<8 | j<< 4 ] = path_to_calib;
        } 


    int succ = 1;
    while( succ ) {
        succ = unpacker::get_time_window( meta_data, original_data, filtered_data, preproc_data, paths_to_tdc_calib, fp );

         /* print original data */
        // for (auto const &pair: original_data)
        // {
        //     printf("{%02x\n", pair.first);

        //     for (auto const &val : pair.second)
        //     {
        //         printf("\t%02x (%d, %d, %.0f),\n", val.sample,
        //                                           val.is_falling_edge,
        //                                           val.channel_id,
        //                                           val.time);
        //     }
            
        //     printf("}\n");
        // }
    
        // /* print filtered data */
        // for (auto const &pair: filtered_data)
        // {
        //     printf("{%02x\n", pair.first);

        //     for (auto const &val : pair.second)
        //     {
        //         printf("\t%02x (%d, %d, %.0f),\n", val.sample,
        //                                           val.is_falling_edge,
        //                                           val.channel_id,
        //                                           val.time);
        //     }
            
        //     printf("}\n");
        // }

        // /* print preproc data */      
        // for (auto const &pair: preproc_data)
        // {
        //     printf("{%02x\n", pair.first);

        //     for (auto const &val : pair.second)
        //     {
        //         printf("\t(%.0f %.0f %d %d),\n", val.lead_time,
        //                                     val.tot_time,
        //                                     val.strip_id,
        //                                     val.multiplicity);
        //     }

        //     printf("}\n");
        // }


        // getchar();
    }


    return 0;
}
