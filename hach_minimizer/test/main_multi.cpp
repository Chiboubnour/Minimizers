#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <cstdint>
#include <string>
#include <inttypes.h> 
#include <cstdlib>
#include <numeric>   
#include <cstring>   
#include <algorithm> 
#include <getopt.h>
#include <omp.h>

#include <filesystem>
namespace fs = std::filesystem;


// Déclaration de la fonction de hachage
extern "C" void minimizer_cpu_multi(
    const uint64_t* input_minimizers,
          uint64_t* output_hashes,
          unsigned int data_size_words);


void read_uint64_t_file(
    const std::string&           filename,
          std::vector<uint64_t>& i_buffer,
    const int nElements )
{
    //
    //
    //
    fs::path filepath(filename);
    if ( fs::exists(filepath) == false )
    {
        std::cout << "(EE) Unable to open file (" << filename << ")" << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //
    //
    //
    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file) {
        std::cout << "(EE) Unable to open file (" << filepath << ")" << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //
    //
    //
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    //
    //
    //
    if (size % sizeof(uint64_t) != 0)
    {
        std::cout << "(EE) The file size is not %8 (uint8_t)" << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }

    //
    //
    //
    int nCases = size / sizeof(uint64_t);
    if( nElements != -1 )
    {
        nCases = (nCases > nElements) ? nElements : nCases;
    }

    
    //
    //
    //
    i_buffer.resize( nCases );
    size_t nBytes = nCases * sizeof(uint64_t);
    if ( !file.read(reinterpret_cast<char*>(i_buffer.data()), nBytes) )
    {
        std::cout << "(EE) An error happens during the data reading" << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    
}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
int main(int argc, char* argv[])
{
    //
    //
    //
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <fichier_minimizers> [num_tests=10]" << std::endl;
        return EXIT_FAILURE;
    }

    int verbose_flag =  0;
    int help_flag    =  0;
    int nElements    = -1;
    int nTests       = 64;
    int nCores       =  1;

    std::string i_file = "";
    std::string o_file = "";

    static struct option long_options[] ={
            {"verbose",     no_argument,       0,  'v'},
            {"help",        no_argument,       0,  'h'},
            {"input",       required_argument, 0,  'i'},
            {"output",      required_argument, 0,  'o'},
            {"cores",       required_argument, 0,  'c'},
            {"nhash",       required_argument, 0,  'n'},
            {"ntest",       required_argument, 0,  't'},
            {0, 0, 0, 0}
    };

    /**
     * 
     * getopt_long stores the option index here.
     * 
     * */
    int c;
    int option_index = 0;
    while( true )
    {
        c = getopt_long (argc, argv, "i:o:n:t:c:vh", long_options, &option_index);
        if (c == -1)
            break;

        switch ( c )
        {
            case 'i': // input file
                i_file = optarg;
                break;

            case 'o': // output file
                o_file = optarg;
                break;

            case 'n':   // #hash
                nElements = 1024 * 1024 * std::atoi(optarg);
                break;

            case 't':   // #tests
                nTests = std::atoi(optarg);
                break;

            case 'c':   // #tests
                nCores = std::atoi(optarg);
                break;

            case 'v':   // verbose
                verbose_flag = true;
                break;

            case 'h':   // help
                help_flag = true;
                break;

            case '?':
                break;

            default:
                abort ();
        }
    }

    //
    //
    //
    if( help_flag == true )
    {
        std::cout << "Usage: " << argv[0] << " -i <minimizer file> -o <ouput file>" << std::endl;
        std::cout << "--input  filename (-i) : " << std::endl;
        std::cout << "--output filename (-o) : " << std::endl;
        std::cout << "--nhash  integer  (-n) : " << std::endl;
        std::cout << "--ntest  integer  (-t) : " << std::endl;
        std::cout << "--help            (-h) : " << std::endl;
        std::cout << "--verbose         (-v) : " << std::endl;
        return EXIT_FAILURE;
    }

    //
    //
    //
    if( i_file.size() == 0 )
    {
        std::cout << "(EE) No input file provided..." << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //
    //
    //
    if( o_file.size() == 0 )
    {
        std::cout << "(EE) No ouput file provided..." << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //
    // 
    //
    if( (nElements < (1024*1024)) && (nElements != -1) )
    {
        std::cout << "(EE) The number of hash to process (" << nElements << ") is (< 1)..." << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //
    // 
    //
    if( nTests < 1 )
    {
        std::cout << "(EE) The number of test to execute is (< 1)..." << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        exit(EXIT_FAILURE);
    }


    //
    // Testing if the file exists, otherwise it is the user's problem
    //
    std::ifstream iTest( i_file );
    if( iTest.is_open() == false )
    {
        std::cout << "(EE) Unable to open file (" << i_file << ")" << std::endl;
        std::cout << "(EE) Issue located (" << __FILE__ << " :: " << __LINE__ << ")" << std::endl;
        return EXIT_FAILURE;
    }
    iTest.close();

    //
    // 
    //
    if( verbose_flag == true )
    {
        std::cout << "#(II) Loading files from file (" << i_file << ")" << std::endl;
    }
    
    std::vector<uint64_t> minimizers;
    read_uint64_t_file(
        i_file,         // input file
        minimizers,     // the buffer for uint64_t values
        nElements );    // the maximum number of elements to load

    //
    //
    //
    std::vector<uint64_t> hashes( minimizers.size() );
    if( verbose_flag == true ){
        std::cout << "#(II) Starting test procedure (" << nTests << " iter.)" << std::endl;
    }

    //
    //
    //
    std::vector<double> durations_us;
    durations_us.reserve(nTests);

    //
    //
    //
    omp_set_num_threads(nCores);
 
    //
    //
    //
    for (int i = 0; i < nTests; ++i)
    {
        //
        // Starting the time measure
        //
        auto start = std::chrono::high_resolution_clock::now();
    
        //
        // We launch the HASH computation
        //
        minimizer_cpu_multi(minimizers.data(), hashes.data(), minimizers.size());
    
        //
        // End of the time measure + data pushback for latter statistics
        //
        auto end       = std::chrono::high_resolution_clock::now();
        double time_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        durations_us.push_back(time_us);

        //
        // We print the N first and N last HASH values for debugging purpose (comparing execution on /= computers)
        // VERBOSE mode only !
        //
        if ( (i == 0) && (verbose_flag == true) )
        {
            std::cout << "#" << std::endl;
            std::cout << "#(II) Premiers hash générés :" << std::endl;
            for (size_t j = 0; j < std::min<size_t>(5, hashes.size()); ++j)
            {
                printf("#(II) - hash[%zu] = %16.16" PRIX64 "\n", j, hashes[j]);

            }
            std::cout << "#" << std::endl;
        }
    
        //
        // A fake computation added here to avoid that GCC/CLANG remove the "useless" procesing
        //
        if (i < nTests - 1) {
            minimizers[rand() % minimizers.size()] = hashes[rand() % hashes.size()];
        }
    
    }
    
    //
    // Execution time statistic computations
    //
    std::sort(durations_us.begin(), durations_us.end());
    double time_min    = durations_us.front() / 1e6;
    double time_max    = durations_us.back() / 1e6;
    double time_median = durations_us[durations_us.size() / 2] / 1e6;
    double time_avg    = std::accumulate(durations_us.begin(), durations_us.end(), 0.0) / (durations_us.size() * 1e6);

    double total_hashes   = static_cast<double>(minimizers.size());
    double hashes_per_sec = total_hashes / time_avg;

    const int    mHash = minimizers.size() / 1e6;
    const int    Mbytes = static_cast<int>(minimizers.size() * sizeof(uint64_t) / (1024.0 * 1024.0));
    const double dHash = (hashes_per_sec / 1e6);
    std::cout << "#(II) Final results" << std::endl;
    std::cout << "#(II) - #of hash       : " <<  minimizers.size()                     / (1024.0 * 1024.0) << " Mhash" << std::endl;
    std::cout << "#(II) - #of bytes      : " << (minimizers.size() * sizeof(uint64_t)) / (1024.0 * 1024.0) << " Mbytes" << std::endl;
    std::cout << "#(II) - #exec. threads : " << nCores   << " core(s)" << std::endl;
    std::cout << "#(II) - Temps moyen    : " << time_avg << " s" << std::endl;
    std::cout << "#(II) - Temps min.     : " << time_min << " s" << std::endl;
    std::cout << "#(II) - Temps max.     : " << time_max << " s" << std::endl;
    std::cout << "#(II) - Temps médian   : " << time_median << " s" << std::endl;
    std::cout << "#(II) - Débit          : " << dHash << " M hash/s" << std::endl;

    printf("%4d %4d  %1.6f %1.6f %1.6f %1.6f %7.1f\n",
        Mbytes, mHash, time_avg, time_min, time_max, time_median, dHash);
        
    return EXIT_SUCCESS;
}
