#pragma once

#include "../files/stream_reader.hpp"
#include "../files/stream_writer.hpp"

class parser_fastq
{
public:
    static bool execute(
        stream_reader&  i_file, 
        stream_writer&  o_file);
};
