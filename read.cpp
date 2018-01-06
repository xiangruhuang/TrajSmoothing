#include<iostream>
#include<fstream>
#include<map>
#include<vector>
#include<assert.h>
#include<algorithm>

using namespace std;

typedef struct _data_unit{
    int CUID; // Car ID (unique in all files)
    int UTC; // UTC time (seconds elapsed since 1970-01-01 00:00 GMT)
    int LAT; // Location Latitude (*1e5)
    int LONG; // Location Longtitude (*1e5)
    int HEAD; // Orientation (0-360)
    int SPEED; // Speed (cm/sec)
    int OCCUPANT; // The Taxi is occupied(1) or empty(0)
} data_unit;

struct Car{
    vector<pair<int, pair<double, double>>> traj;
    int car_id;
};

data_unit read_one_data_unit(ifstream* fin){
    _data_unit dat;
    
    fin->read((char*)&dat.CUID, sizeof(dat.CUID));
    fin->read((char*)&dat.UTC, sizeof(dat.UTC));
    fin->read((char*)&dat.LAT, sizeof(dat.LAT));
    fin->read((char*)&dat.LONG, sizeof(dat.LONG));
    fin->read((char*)&dat.HEAD, sizeof(dat.HEAD));
    fin->read((char*)&dat.SPEED, sizeof(dat.SPEED));
    fin->read((char*)&dat.OCCUPANT, sizeof(dat.OCCUPANT));
    //cout << "CUID=" << dat.CUID << ", UTC=" << dat.UTC << ", LAT=" << dat.LAT;
    //cout << ", LONG=" << dat.LONG << ", HEAD=" << dat.HEAD << ", SPEED=" << dat.SPEED;
    //cout << ", OCCUPANT=" << dat.OCCUPANT;
    //cout << endl;
    
    return dat;
}

// Return Dense Unique ID of this car, register to carid_map if necessary
//
int register_car_or_not(data_unit dat, map<int, int>& carid_map, vector<Car>& car){ 
    map<int, int>::iterator it = carid_map.find(dat.CUID);
    int idx; // dense unique id of this car
    if (it == carid_map.end()){
        idx = car.size();
        carid_map.insert(it, make_pair(dat.CUID, idx));
        Car new_car;
        new_car.car_id = dat.CUID;
        car.push_back(new_car);
    } else {
        idx = it->second;
    }
    return idx;
}

int main(int argc, char** agrv){


    vector<Car> car;

    map<int, int> carid_map;
    vector<pair<double, pair<double, double>>> traj; // <time, <LAT, LONG>>

    ofstream fout;
    fout.open("real/trajectories.txt", fstream::out);
    for (int i = 1; i <= 29; i++){
        ifstream fin;
        string prefix = "data/output_05";
        string suffix = ".dat";
        if (i == 10 || i == 20){
            //no data for these two days
            continue;
        }
        string date = to_string(i);
        if (i < 10){
            date = "0" + date;
        }
        string filename = prefix+date+suffix;
        cout << "reading from " << filename << endl;
        fin.open(filename, fstream::in | fstream::binary);

        int num_sample = 0;
        while (!fin.eof()){
            data_unit dat = read_one_data_unit(&fin);
            int idx = register_car_or_not(dat, carid_map, car);
            car[idx].traj.push_back(make_pair(dat.UTC, make_pair(dat.LAT, dat.LONG)));
            //num_sample++;
            //cout << "counter=" << counter << ", #sample=" << num_sample << endl;
        }
        for (int c = 0; c < car.size(); c++){
            Car& car_c = car[c];
            sort(car_c.traj.begin(), car_c.traj.end(), less<pair<int, pair<double, double>>>());
            int last_time = -1000000;
            for (int k = 0; k < car_c.traj.size(); k++){
                if (k != 0){
                    fout << " ";
                }
                int UTC = car_c.traj[k].first;
                double LAT = car_c.traj[k].second.first/1e5;
                double LONG = car_c.traj[k].second.second/1e5;
                fout << LONG << " " << LAT;
                //cout << "last_time=" << last_time << ", UTC=" << UTC << endl;
                assert(last_time <= UTC);
                last_time = UTC;
            }
            fout << endl;
            car_c.traj.clear();
        }
        fin.close();
    }
    fout.close();


}
