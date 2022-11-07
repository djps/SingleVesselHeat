/**
 @file heatReader.h
*/

#ifndef HEATREADER_H
#define HEATREADER_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cerrno>

#include <iomanip>

#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <exception>

#include "rapidxml.hpp"
#include "types.h"

#include "json.hpp"

using namespace std;
using namespace rapidxml;
using nlohmann::json;

static inline std::string &ltrim(std::string &s);

static inline std::string &rtrim(std::string &s);

static inline std::string &trim(std::string &s);

HeatInputData traverse_xml(const string &input_xml);

HeatInputData traverse_json(const string &input_json);

#endif
