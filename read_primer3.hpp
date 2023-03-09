#ifndef READPRIMER3_H
#define READPRIMER3_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "libprimervg.hpp"

int read_primer_entry(std::ifstream* file, vgprimers* primers);

#endif