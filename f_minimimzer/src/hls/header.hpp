#pragma once
//
#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>
//
//  Parametres globaux du noyau. PAR_FACTOR n'est PLUS un global : c'est un
//  parametre TEMPLATE de chaque module (thr_xxx<PAR_FACTOR>), pour generer des
//  circuits a differents facteurs de parallelisme (8, 16, 32, 64...).
//
constexpr int SMER_SIZE   = 19;   // taille d'un s-mer (bases)
constexpr int WINDOW_SIZE = 16;   // taille de la fenetre du minimizer (s-mers)
//
