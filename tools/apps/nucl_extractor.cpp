#include <cstdio>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <sys/stat.h>
#include <dirent.h>
//
//
#include "../src/files/stream_reader.hpp"
#include "../src/files/stream_reader_library.hpp"
#include "../src/files/stream_writer_library.hpp"
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
uint64_t get_file_size(const std::string& filen) {
    struct stat file_status;
    if (stat(filen.c_str(), &file_status) < 0) {
        return -1;
    }
    return file_status.st_size;
}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
int main(int argc, char *argv[])
{
    std::string i_file;
    std::string o_file;

    int nElements    = 0;
    int verbose_flag = 0;
    int help_flag    = 0;

    static struct option long_options[] ={
            {"verbose",     no_argument, &verbose_flag, 1},
            {"help",        no_argument, 0, 'h'},
            {"ifile",       required_argument, 0,  'i'},
            {"ofile",       required_argument, 0,  'o'},
            {"size",        required_argument, 0,  'n'},
            {0, 0, 0, 0}
    };

 
    int option_index = 0;
    int c;
    while( true )
    {
        c = getopt_long (argc, argv, "i:o:vh", long_options, &option_index);
        if (c == -1)
            break;

        switch ( c )
        {
            case 'i':
                i_file = optarg;
                break;

            case 'o':
                o_file = optarg;
                break;

            case 's':
                nElements = std::atoi(optarg);
                break;

            case 'v':
                verbose_flag = true;
                break;

            case 'h':
                help_flag = true;
                break;

            case '?':
                break;

            default:
                abort ();
        }
    }

   
    if( i_file.size() == 0 )
        printf ("(EE) First file to merge is undefined !\n");

    if ( (optind < argc) || (help_flag == true) || (i_file.size() == 0) )
    {
        printf ("(II) Usage :\n");
        printf ("(II) ./merge -f <first file> -s <second file> -o <output file> -c <first color> -C <second color>");
        printf ("(II)\n");
        printf ("(II) Options :\n");
        printf ("(II) - \n");
        putchar ('\n');
        exit( EXIT_FAILURE );
    }


    const int64_t i_size_bytes  = get_file_size(i_file);
    const int64_t i_size_Kbytes = (uint64_t)i_size_bytes / 1024;
    const int64_t i_size_Mbytes = i_size_Kbytes / 1024;
    if( i_size_bytes == -1 )
    {
        printf("(EE) The first file to merge does not exist (%s)\n", i_file.c_str());
        printf("(EE) Error was detected in %s (line %d)\n", __FILE__, __LINE__);
        exit( EXIT_FAILURE );
    }

 
    printf("(II)\n");
    printf("(II) Input document information\n");
    printf("(II) - filename    : %s\n", i_file.c_str());
    printf("(II) - filesize    : %llu bytes\n", i_size_bytes);
    if( i_size_Kbytes <= 32 )
        printf("(II) -             : %llu Kbytes\n", i_size_Kbytes);
    else
        printf("(II) -             : %llu Mbytes\n", i_size_Mbytes);
    printf("(II)\n");

//    double start_time = omp_get_wtime();

    stream_reader* reader = stream_reader_library::allocate( i_file );
    stream_writer* writer = stream_writer_library::allocate( o_file );

    // virtual bool is_open () = 0;
    // virtual bool is_open () = 0;

    uint8_t* i_buffer = new uint8_t[65536];
    uint8_t* o_buffer = new uint8_t[16384];
    
//    virtual void close   () = 0;
//    virtual bool is_eof  () = 0;
//    virtual int  read    (void* buffer, int eSize, int eCount) = 0;

    while (true)
    {
        const int nRead = reader->read(i_buffer, sizeof(uint8_t), 65536 );

        //
        // Si il n'y a plus de données dans le fichier, on s'arrete !
        //
        if( reader->is_eof() == true ){
            break;
        }
    }

    delete[] i_buffer;
    delete[] o_buffer;
    
    reader->close();
    writer->close();
    
    delete reader;
    delete writer;

    return 0;
}