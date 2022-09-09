#include <vector>

class BlackScholesModel
{
private:
    float strike;
    string optionType;
    vector<double> payoffCoeff;
    int timestepNumber;
    int hedgingDates;
    float fdStep;
    long sampleNumber;

public:
    BlackScholesModel(/* args */);
    ~BlackScholesModel();
};

BlackScholesModel::BlackScholesModel(/* args */)
{
}

BlackScholesModel::~BlackScholesModel()
{
}
