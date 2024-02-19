#pragma once
#ifndef _CONSOLE_COLOURS_H
#define _CONSOLE_COLOURS_H

#include <string>

namespace Console {
	std::string Reset = "\x1b[0m";
	std::string Bright = "\x1b[1m";
	std::string Dim = "\x1b[2m";
	std::string Underscore = "\x1b[4m";
	std::string Blink = "\x1b[5m";
	std::string Reverse = "\x1b[7m";
	std::string Hidden = "\x1b[8m";

	namespace Foreground {
		std::string Black = "\x1b[30m";
		std::string Red = "\x1b[31m";
		std::string Green = "\x1b[32m";
		std::string Yellow = "\x1b[33m";
		std::string Blue = "\x1b[34m";
		std::string Magenta = "\x1b[35m";
		std::string Cyan = "\x1b[36m";
		std::string White = "\x1b[37m";
	}

	namespace Background {
		std::string Black = "\x1b[40m";
		std::string Red = "\x1b[41m";
		std::string Green = "\x1b[42m";
		std::string Yellow = "\x1b[43m";
		std::string Blue = "\x1b[44m";
		std::string Magenta = "\x1b[45m";
		std::string Cyan = "\x1b[46m";
		std::string White = "\x1b[47m";
	}

	using namespace Foreground;
}

#endif