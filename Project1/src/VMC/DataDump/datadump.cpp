#include "datadump.h"

/*
This class is made to make it easy to dump data. By now this class uses a template
which will accept either double, int or a vector of doubles.

This class is used by simply making an instance which take the name of the file
you want to contain the data. The constructor can also take the filename of a stamp file.
A stamp is an additional variable one can save together with the main data (optional and
somewhat pointless).

Use:

DataDump<type> instance_name("filename");

instance_name.push_back(data_point); //Push back a new data point
instance_name.push_back_stamp(data_stamp); //Push back a new stamp

instance_name.dump(data_point); //Directly dumps a data point to file(slow and should not be used)

instance_name.dump_all(data_point); //Dumps everthing, should be used at the end

instance_name.dump_metadata(); //Dumps all the parameters from parameter file as meta data

*/

template<class T>
void DataDump<T>::push_back(T data_point){
    data.push_back(data_point);
}

template<class T>
void DataDump<T>::push_back_stamp(double data_point){
    data_stamp.push_back(data_point);
}

template<class T>
void DataDump<T>::dump(T data_point){
    outfile << data_point;
}

template<>
void DataDump<std::vector<double>>::dump(std::vector<double> data_point){
    dump_vector(data_point);
}

template<class T>
void DataDump<T>::dump(T data_point,double stamp){
    dump(data_point);
    if(include_stamp){
        stampfile << stamp;
    }

}


template<>
void DataDump<std::vector<double>>::dump(std::vector<double> data_point,double stamp){
    dump_vector(data_point);
    if(include_stamp){
        stampfile << stamp;
    }

}

template<class T>
void DataDump<T>::dump_all(){
    int data_size = data.size();
    for(int i = 0; i<data_size;i++){
        outfile << data[i] << " ";
    }
    if(include_stamp){
        int stamp_size = data_stamp.size();
        for(int i = 0; i<   stamp_size;i++){
            stampfile << data_stamp[i] << " ";
        }
    }
}

template<>
void DataDump<std::vector<double>>::dump_all(){
    int data_size = data.size();
    for(int i = 0; i<data_size;i++){
        dump_vector(data[i]);
    }
    if(include_stamp){
        int stamp_size = data_stamp.size();
        for(int i = 0; i<   stamp_size;i++){
            stampfile << data_stamp[i] << " ";
        }
    }
}



template<class T>
void DataDump<T>::dump_vector(std::vector<double> data){
    int size = data.size();
    for(int i = 0; i < size; i++){
        outfile << data[i] << "\n";
    }
}

template<class T>
void DataDump<T>::dump_metadata(std::string m_location){
    std::fstream metafile(m_location,std::fstream::out);
    metafile << "MC_cycle " << Parameters::MC_cycles << "\n";

    metafile << "Alpha_min " << Parameters::alpha_min << "\n";
    metafile << "Alpha_max " << Parameters::alpha_max << "\n";
    metafile << "Alpha_num " << Parameters::alpha_num << "\n";

    metafile << "Beta " << Parameters::beta << "\n";
    metafile << "Omega " << Parameters::omega << "\n";
    metafile << "Omega_z " << Parameters::omega_z << "\n";


    metafile << "Dimensions " << Parameters::dimension << "\n";
    metafile << "N " << Parameters::N << "\n";

    metafile << "a " << Parameters::a << "\n";
    metafile << "D " << Parameters::D << "\n";
    metafile << "dx " << Parameters::dx << "\n";

    metafile.close();

}

template class DataDump<double>;
template class DataDump<int>;
template class DataDump<std::string>;
template class DataDump<std::vector<double>>;





