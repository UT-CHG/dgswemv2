#ifndef NUMA_CONFIGURATION_HPP
#define NUMA_CONFIGURATION_HPP

class NumaConfiguration
{
public:
    NumaConfiguration()= default;

    NumaConfiguration(const std::string& type)
    {
        if ( type.compare("default") ) {
            numa_wghts = { 1. };
        } else if ( type.compare("stampede") ) {
            numa_wghts = { 1. , 1. };
        } else {
            std::cerr << "Numa Configuration: " << type
                      << " is invalid\n";
            std::cerr << "Aborting\n";
            exit(1);
        }
    }

    inline
    const std::vector<double>& get_numa_wghts() const
    { return numa_wghts; }

    inline
    std::size_t get_num_numa_domains() const
    { return numa_wghts.size(); }

private:
    std::vector<double> numa_wghts;
};
#endif
