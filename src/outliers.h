// calculate the mean and standard deviation
void cal_mean_std(std::vector<double> &value, double &mean, double &std)
{
    double sum = 0;
    double variance = 0;
    assert(value.size() > 1);
    for (int i = 0; i < value.size(); i++)
    {
        sum += value[i];
    }
    mean = sum / value.size();
    for (int i = 0; i < value.size(); i++)
    {
        variance += pow((value[i] - mean), 2);
    }
    variance = variance / (value.size() - 1);
    std = sqrt(variance);
}

void remove_outliers(std::vector<double> &with_outliers,
                     std::vector<double> &without_outliers)
{
    double mean;
    double std;
    cal_mean_std(with_outliers, mean, std);
    double upper_limit = mean + std;
    double lower_limit = mean - std;

    for (int i = 0; i < with_outliers.size(); i++)
    {
        if (with_outliers[i] <= upper_limit && with_outliers[i] >= lower_limit)
        {
            without_outliers.push_back(with_outliers[i]);
        }
    }
}
