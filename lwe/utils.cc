#include "lwe/utils.h"

int small_modexp_2n(int base, int exp, int mod) {
  int result = 1;
  for (int i = 0; i < exp; i++) {
    result = (result * base) % (2 * mod); // 2 * mod because in X^n + 1 cylotomic polynomial, we have that X^{2*n} = 1 unlike X^n = -1.
  }
  return result;
}

std::pair<int, int> find_optimal_bsgs(int rotations, int blocks) {
  int optimal_cost = INT_MAX;
  int optimal_giant = 0;
  int optimal_baby = 0;
  for (int giant = 1; giant <= rotations; giant++) {
    int baby = std::ceil(((float)rotations) / giant);
    int cost = baby + (giant * blocks);
    if (cost < optimal_cost) {
      optimal_cost = cost;
      optimal_baby = baby;
      optimal_giant = giant;
    }
  }
  //std::cout << "Rotations: " << rotations << " Blocks: " << blocks << " Optimal Baby: " << optimal_baby << " Optimal Giant: " << optimal_giant << " Optimal cost: " << optimal_cost << std::endl;
  return std::make_pair(optimal_baby, optimal_giant);
}

int DivideAndRoundUp(int x, int y) {
    return (x + y - 1) / y; 
}