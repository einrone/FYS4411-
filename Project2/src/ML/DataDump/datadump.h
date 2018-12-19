#ifndef DATADUMP_H
#define DATADUMP_H
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../Parameters/parameters.h"

template<class T>
class DataDump
{
public:
    DataDump(std::string m_location, std::string stamp_location)
        :DataDump(m_location)
    {
        stampfile.open(stamp_location,std::fstream::out);
        include_stamp = true;
    }

    DataDump(std::string m_location){
        outfile.open(m_location,std::fstream::out | std::fstream::binary);
    }
    ~DataDump(){
        outfile.close();
        if(include_stamp)
            stampfile.close();
    }

    void push_back(T);

    void dump(T);
    void dump(T,double);

    void dump_vector(std::vector<double>);

    void dump_all();

    void dump_metadata(std::string m_location);


    std::vector<T> data;
    std::vector<double> data_stamp;

    void push_back_stamp(double data_point);
private:
    std::fstream outfile;
    std::fstream stampfile;

    bool include_stamp;

};

#endif // DATADUMP_H
