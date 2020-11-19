#include "CLI11.hpp"


#ifndef KFF_TOOLS
#define KFF_TOOLS


class KffTool {
public:
	virtual void cli_prepare(CLI::App * subapp) = 0;
};


#endif
