#include "CLI11.hpp"
#include "kff-cpp-api/kff_io.hpp"


#ifndef KFF_TOOLS
#define KFF_TOOLS


class KffTool {
public:
	CLI::App * subapp = nullptr;
	virtual void cli_prepare(CLI::App * subapp) = 0;
	virtual void exec() = 0;
};


#endif
