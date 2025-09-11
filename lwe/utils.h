#include <iostream>
#include <utility>
#include <cmath>
#include <climits>

int small_modexp_2n(int base, int exp, int mod);

std::pair<int, int> find_optimal_bsgs(int rotations, int blocks);

int DivideAndRoundUp(int x, int y);