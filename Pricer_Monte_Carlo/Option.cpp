#include <vector>

class Option
{
private:
    int optionSize;
    vector<double> spot;
    float maturity;
    vector<float> volatility;
    float interestRate;
    float correlation;
    vector<double> trend;

public:
    Option(/* args */);
    ~Option();
};

Option::Option(/* args */)
{
}

Option::~Option()
{
}
