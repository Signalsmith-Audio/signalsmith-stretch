#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <cstdlib> // exit() and codes

#include "console-colours.h"

/** Expected use:

		SimpleArgs args(argc, argv);

		// positional argument
		std::string foo = args.arg<std::string>("foo");
		// optional argument
		std::string bar = args.arg<std::string>("bar", "a string for Bar", "default");
		// --flag=value
		double = args.flag<double>("baz", "an optional flag", 5);
		
		// Exits if "foo" not supplied
		args.errorExit();

	If you have multiple commands, each with their own options:

		// Switches based on a command
		if (args.command("bink", "Bink description")) {
			// collect arguments for the command
		}
		// Exits with a help message (and list of commands) if no command matched
		args.errorCommand();
		
	By default, a flag of "-h" (or a command of "help", if any commands are used) prints a help message.  To override:
		SimpleArgs args(argc, argv);
		args.helpFlag("h");
		args.helpCommand("help");
	
**/
class SimpleArgs {
	int argc;
	const char* const* argv;

	template<typename T>
	T valueFromString(const char *arg);
	
	std::string parsedCommand;
	struct Keywords {
		std::string keyword;
		std::string description;
		bool isHelp;
	};
	std::vector<Keywords> keywordOptions;
	std::vector<Keywords> argDetails;
	std::vector<Keywords> flagOptions;
	std::set<std::string> flagSet;
	void clearKeywords() {
		keywordOptions.resize(0);
		flagSet.clear();
		flagOptions.clear();
	}
	
	bool helpMode = false;
	bool checkedHelpCommand = false;
	bool hasError = false;
	std::string errorMessage;
	void setError(std::string message) {
		if (!hasError) {
			hasError = true;
			errorMessage = message;
		}
	}

	std::map<std::string, std::string> flagMap;
	void consumeFlags() {
		while (index < argc && std::strlen(argv[index]) > 1 && argv[index][0] == '-') {
			const char* arg = argv[index++];
			size_t length = strlen(arg);
			
			size_t keyStart = 1, keyEnd = keyStart + 1;
			size_t valueStart = keyEnd;
			// If it's "--long-arg" format
			if (length > 1 && arg[1] == '-') {
				keyStart++;
				while (keyEnd < length && arg[keyEnd] != '=') {
					keyEnd++;
				}
				valueStart = keyEnd;
				if (keyEnd < length) valueStart++;
			}

			std::string key = std::string(arg + keyStart, keyEnd - keyStart);
			std::string value = std::string(arg + valueStart);
			
			if (key == "help") {
				helpMode = true;
			}

			flagMap[key] = value;
	 	}
	 }
	
	int index = 1;
public:
	SimpleArgs(int argc, const char* const argv[]) : argc(argc), argv(argv) {
		std::string cmd = argv[0];
		size_t slashPos = cmd.find_last_of("\\/");
		if (slashPos != std::string::npos) cmd = cmd.substr(slashPos + 1);
		parsedCommand = cmd;
	}

	void help(std::ostream& out=std::cerr) const {
		std::string parsedCommand = this->parsedCommand;
		if (keywordOptions.size() > 0) {
			parsedCommand += std::string(" <command>");
		}
		out << "Usage:\n\t" <<  parsedCommand << "\n\n";
		if (keywordOptions.size() > 0) {
			out << "Commands:\n";
			for (unsigned int i = 0; i < keywordOptions.size(); i++) {
				out << "\t" << keywordOptions[i].keyword;
				if (keywordOptions[i].isHelp) out << " [command...]";
				if (keywordOptions[i].description.size()) out << "  -  " << keywordOptions[i].description;
				out << "\n";
			}
			out << "\n";
		}
		if (argDetails.size() > 0) {
			out << "Arguments:\n";
			for (Keywords const &arg : argDetails) {
				out << "\t" << arg.keyword;
				if (arg.description.size()) out << "  -  " << arg.description;
				out << "\n";
			}
			out << "\n";
		}
		if (flagOptions.size() > 0) {
			out << "Options: " << Console::Dim << "(--arg=value)" << Console::Reset << "\n";
			for (Keywords const &pair : flagOptions) {
				out << "\t" << (pair.keyword.length() > 1 ? "--" : "-") << pair.keyword;
				if (pair.description.size()) out << "  -  " << pair.description;
				out << "\n";
			}
			out << "\n";
		}
	}
	
	bool isHelp() const {
		return helpMode;
	}
	bool finished() const {
		return index >= argc;
	}
	std::string peek() const {
		return (index >= argc) ? "" : argv[index];
	}

	int errorExit(std::ostream& out=std::cerr) const {
		if (hasError || helpMode) {
			help(out);
			if (!helpMode) {
				out << Console::Red << errorMessage << Console::Reset << "\n";
			}
			std::exit((!helpMode && hasError) ? EXIT_FAILURE : EXIT_SUCCESS);
		}
		return 0;
	}
	int errorExit(std::string forcedError, std::ostream& out=std::cerr) const {
		if (hasError) return errorExit(out); // Argument errors take priority
		out << Console::Red << forcedError << Console::Reset << "\n";
		std::exit(EXIT_FAILURE);
		return 0;
	}
	int errorCommand(std::string message="", std::ostream& out=std::cerr) const {
		if (keywordOptions.size()) {
			// We expected a command, but didn't match on any
			if (helpMode) return errorExit(out);
			if (index >= argc) help(out);
			if (message.length() == 0) {
				message = (index < argc) ? std::string("Unknown command: ") + argv[index] : "Missing command";
			}
			errorExit(message, out);
		}
		return 0;
	}
	
	template<typename T=std::string>
	T arg(std::string name, std::string longName, T defaultValue) {
		consumeFlags();
		if (index < argc) clearKeywords();
		parsedCommand += std::string(" [") + name + "]";
		argDetails.push_back(Keywords{name, longName, false});

		if (index >= argc) return defaultValue;
		return valueFromString<T>(argv[index++]);
	}

	template<typename T=std::string>
	T arg(std::string name, std::string longName="") {
		consumeFlags();
		if (index < argc) clearKeywords();
		parsedCommand += std::string(" <") + name + ">";
		argDetails.push_back(Keywords{name, longName, false});

		if (index >= argc) {
			if (longName.length() > 0) {
				setError("Missing " + longName + " <" + name + ">");
			} else {
				setError("Missing argument <" + name + ">");
			}
			return T();
		}

		return valueFromString<T>(argv[index++]);
	}

	bool command(std::string keyword, std::string description="", bool isHelp=false) {
		consumeFlags();
		if (index == 1) {
			helpCommand();
		}
		if (index < argc && !keyword.compare(argv[index])) {
			clearKeywords();
			index++;
			if (!isHelp) parsedCommand += std::string(" ") + keyword;
			return true;
		}
		keywordOptions.push_back(Keywords{keyword, description, isHelp});
		return false;
	}
	bool helpCommand(std::string keyword="help") {
		if (!checkedHelpCommand && index == 1) {
			keywordOptions.push_back(Keywords{keyword, "", true});
			if (index < argc && !keyword.compare(argv[index])) {
				index++;
				helpMode = true;
			}
		}
		checkedHelpCommand = true;
		return helpMode;
	}

	template<typename T=std::string>
	T flag(std::string key, std::string description, T defaultValue) {
		consumeFlags();
		if (!hasFlag(key, description)) return defaultValue;

		auto iterator = flagMap.find(key);
		return valueFromString<T>(iterator->second.c_str());
	}
	template<typename T=std::string>
	T flag(std::string key, T defaultValue) {
		consumeFlags();
		if (!hasFlag(key, "")) return defaultValue;

		auto iterator = flagMap.find(key);
		return valueFromString<T>(iterator->second.c_str());
	}
	template<typename T=std::string>
	T flag(std::string key) {
		return flag<T>(key, T());
	}
	bool hasFlag(std::string key, std::string description="") {
		consumeFlags();
		auto iterator = flagSet.find(key);
		if (iterator == flagSet.end()) {
			flagSet.insert(key);
			flagOptions.push_back(Keywords{key, description, false});
		} else if (description.length() > 0) {
			bool found = false;
			for (auto &option : flagOptions) {
				if (option.keyword == key) {
					option.description = description;
					found = true;
					break;
				}
			}
			if (!found) {
				flagOptions.push_back(Keywords{key, description, false});
			}
		}

		auto mapIterator = flagMap.find(key);
		return mapIterator != flagMap.end();
	}
	bool helpFlag(std::string key, std::string description="shows this help") {
		consumeFlags();
		hasFlag(key, description);
		auto iterator = flagMap.find(key);
		helpMode = (iterator != flagMap.end());
		return helpMode;
	}
};

template<>
std::string SimpleArgs::valueFromString(const char *arg) {
	return arg;
}
template<>
const char * SimpleArgs::valueFromString(const char *arg) {
	return arg;
}
template<>
int SimpleArgs::valueFromString(const char *arg) {
	return std::stoi(arg);
}
template<>
long SimpleArgs::valueFromString(const char *arg) {
	return std::stol(arg);
}
template<>
unsigned long SimpleArgs::valueFromString(const char *arg) {
	return std::stoul(arg);
}
template<>
float SimpleArgs::valueFromString(const char *arg) {
	return std::stof(arg);
}
template<>
double SimpleArgs::valueFromString(const char *arg) {
	return std::stod(arg);
}
