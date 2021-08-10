// https://cplusplus.com/reference/random/mersenne_twister_engine/operator()/
// mersenne_twister_engine::operator()
#include <iostream>
#include <chrono>
#include <random>

int main ()
{
  // obtain a seed from the system clock:
  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = std::random_device{}(); // https://stackoverflow.com/questions/24334012/best-way-to-seed-mt19937-64-for-monte-carlo-simulations

  seed *= std::chrono::system_clock::now().time_since_epoch().count();

  std::mt19937 generator (seed);  // mt19937 is a standard mersenne_twister_engine
  for (int i = 0; i < 15; ++i)
    std::cout << "Random value: " << generator()/double(generator.max()) << std::endl;


  std::cout << generator.min() << " and " << generator.max();
  return 0;
}