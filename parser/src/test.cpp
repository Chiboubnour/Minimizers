/*
 *  Copyright (c) 2022 Bertrand LE GAL
 *
 *  This software is provided 'as-is', without any express or
 *  implied warranty. In no event will the authors be held
 *  liable for any damages arising from the use of this software.
 *
 *  Permission is granted to anyone to use this software for any purpose,
 *  including commercial applications, and to alter it and redistribute
 *  it freely, subject to the following restrictions:
 *
 *  1. The origin of this software must not be misrepresented;
 *  you must not claim that you wrote the original software.
 *  If you use this software in a product, an acknowledgment
 *  in the product documentation would be appreciated but
 *  is not required.
 *
 *  2. Altered source versions must be plainly marked as such,
 *  and must not be misrepresented as being the original software.
 *
 *  3. This notice may not be removed or altered from any
 *  source distribution.
 *
 */
#include <cstdio>
#include <cstdlib>
#include <stdbool.h>
#include <cmath>
#include <sys/time.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
//
//
#include "./catch2v3/catch_amalgamated.hpp"
//
//
#include "./tools/compare.hpp"
//
//
#include "./hls/fastq_parser.hpp"
//
#include "./hls/fasta_parser.hpp"
#include "./hls/fasta_parser_neon.hpp"
#include "./hls/fasta_parser_neon_ultimate.hpp"
#include "./hls/fasta_parser_neon_json.hpp"
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
TEST_CASE( "NAIVE", "[naive-test-1]" )
{
    const std::string file_i = "../data/fasta/test_1/input.fasta";
    const std::string file_r = "../data/fasta/test_1/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[naive-test-2]" )
{
    const std::string file_i = "../data/fasta/test_2/input.fasta";
    const std::string file_r = "../data/fasta/test_2/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[naive-test-3]" )
{
    const std::string file_i = "../data/fasta/test_3/input.fasta";
    const std::string file_r = "../data/fasta/test_3/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[naive-test-4]" )
{
    const std::string file_i = "../data/fasta/test_4/input.fasta";
    const std::string file_r = "../data/fasta/test_4/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
TEST_CASE( "NAIVE", "[neon-test-1]" )
{
    const std::string file_i = "../data/fasta/test_1/input.fasta";
    const std::string file_r = "../data/fasta/test_1/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[neon-test-2]" )
{
    const std::string file_i = "../data/fasta/test_2/input.fasta";
    const std::string file_r = "../data/fasta/test_2/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[neon-test-3]" )
{
    const std::string file_i = "../data/fasta/test_3/input.fasta";
    const std::string file_r = "../data/fasta/test_3/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[neon-test-4]" )
{
    const std::string file_i = "../data/fasta/test_4/input.fasta";
    const std::string file_r = "../data/fasta/test_4/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
TEST_CASE( "NAIVE", "[ultra-test-1]" )
{
    const std::string file_i = "../data/fasta/test_1/input.fasta";
    const std::string file_r = "../data/fasta/test_1/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_ultimate_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[ultra-test-2]" )
{
    const std::string file_i = "../data/fasta/test_2/input.fasta";
    const std::string file_r = "../data/fasta/test_2/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_ultimate_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[ultra-test-3]" )
{
    const std::string file_i = "../data/fasta/test_3/input.fasta";
    const std::string file_r = "../data/fasta/test_3/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_ultimate_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[ultra-test-4]" )
{
    const std::string file_i = "../data/fasta/test_4/input.fasta";
    const std::string file_r = "../data/fasta/test_4/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_ultimate_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
TEST_CASE( "NAIVE", "[json-test-1]" )
{
    const std::string file_i = "../data/fasta/test_1/input.fasta";
    const std::string file_r = "../data/fasta/test_1/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_json_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[json-test-2]" )
{
    const std::string file_i = "../data/fasta/test_2/input.fasta";
    const std::string file_r = "../data/fasta/test_2/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_json_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[json-test-3]" )
{
    const std::string file_i = "../data/fasta/test_3/input.fasta";
    const std::string file_r = "../data/fasta/test_3/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_json_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[json-test-4]" )
{
    const std::string file_i = "../data/fasta/test_4/input.fasta";
    const std::string file_r = "../data/fasta/test_4/golden.txt";
    const std::string file_o = "./resu.txt";

    fasta_neon_json_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//
TEST_CASE( "NAIVE", "[naiveq-test-1]" )
{
    const std::string file_i = "../data/fastq/test_1/input.fastq";
    const std::string file_r = "../data/fastq/test_1/golden.txt";
    const std::string file_o = "./resu.txt";

    fastq_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[naiveq-test-2]" )
{
    const std::string file_i = "../data/fastq/test_2/input.fastq";
    const std::string file_r = "../data/fastq/test_2/golden.txt";
    const std::string file_o = "./resu.txt";

    fastq_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[naiveq-test-3]" )
{
    const std::string file_i = "../data/fastq/test_3/input.fastq";
    const std::string file_r = "../data/fastq/test_3/golden.txt";
    const std::string file_o = "./resu.txt";

    fastq_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
TEST_CASE( "NAIVE", "[naiveq-test-4]" )
{
    const std::string file_i = "../data/fastq/test_4/input.fastq";
    const std::string file_r = "../data/fastq/test_4/golden.txt";
    const std::string file_o = "./resu.txt";

    fastq_branchless_parser(file_i, file_o);

    const bool isOK = fichiersEquivalents(file_o, file_r);

    REQUIRE( isOK == true );
}
//
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
//