#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>
#include <utility>

inline
int squashIndex(int i, const std::vector<std::pair<int, int>> & ranges)
{
	int squashed = i;
	for (auto range : ranges) {
		if (i > range.second)
			squashed -= range.second - range.first;
		else if ((i <= range.second) && (i >= range.first))
			return squashed - (i - range.first);
	}
	return squashed;
}

#endif
