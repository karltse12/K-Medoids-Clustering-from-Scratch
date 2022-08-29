#include <iostream>
#include <sstream>
#include <fstream>	
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
#include <cmath> 
#include <algorithm>
using namespace std;
using std::vector; using std::fixed;
using std::setprecision;

// Create variables to store data
	string src_addr;
    string src_port;
    string dst_addr;
    string dst_port;
    string protocol;
    double arrival_time;
    double packet_length;
    string label_name;	

    vector<string> src_addr_vec;
    vector<string> src_port_vec;
    vector<string> dst_addr_vec;
    vector<string> dst_port_vec;
    vector<string> protocol_vec;
    vector<double> arrival_time_vec;
    vector<double> packet_length_vec;
    

// Return vector of row number for a specific cluster
vector<int> cluster_int (string str, string delimiter = " ") {
    vector<int> row_num_vec;
    int start = 0;  
    int end = str.find(delimiter);
    int i;
    int j;
    while (end != -1) {
        std::istringstream(str.substr(start, end - start)) >> i;
        row_num_vec.push_back(i);
        start = end + delimiter.size();
        end = str.find(delimiter, start);
    }
    std::istringstream(str.substr(start, end - start)) >> j;
    row_num_vec.push_back(j);
    return row_num_vec;
}

// Function to calculate the avergae packet length within specific cluster
double avg_length(string flow_row) {

    vector<int> flow = cluster_int(flow_row);
    double total_length =0.0;
    
    for(int i=0; i< flow.size(); i++){
        total_length += packet_length_vec[flow[i]]*100.0/100.0;
    }
    return (total_length/flow.size());
}

// Function to calculate the avergae transferring time within specific cluster
double avg_time(string flow_row) {

    vector<int> flow = cluster_int(flow_row);
    double total_time =0.0;
    
    for(int i= flow.size()-2; i > 0; i--){
        total_time += (arrival_time_vec[flow[i]]-arrival_time_vec[flow[i-1]])*100.0/100.0;
    }
    
    total_time += (arrival_time_vec[flow[0]]-arrival_time_vec[flow[flow.size()-1]])*100.0/100.0;

    return (total_time/(flow.size()-1));
}

// Function to check if the number exist in the vector
bool contain(vector<int> vec, int number) {
    bool result = false;
    if( find(vec.begin(), vec.end(), number) != vec.end()) {
        result = true;
    }
    return result;
}

// Functions to calculate absolute error for pair of data points by Manhattan Distance
double get_AbsoluteError(vector<string> flow_all_str, vector<int> initial_medoids_vec){
    double absolute_error=0;
    vector<int>closest_medoids_vec;
    for(int i=0; i<flow_all_str.size();i++) {
    // excldue initial medoids itselves
    if(contain(initial_medoids_vec, i)){
        // -1 means no medoids will be assigned since it is an initial mediod already
        closest_medoids_vec.push_back(i);  
    }
        else {
            double distance = 1000000;     // Set an extreme high default value for distance
            int cloest_medoids;
            for(int j=0; j<initial_medoids_vec.size(); j++) {
                    double cal_distance = abs((avg_time(flow_all_str[i])-avg_time(flow_all_str[initial_medoids_vec[j]]))) 
                                        + abs((avg_length(flow_all_str[i])-avg_length(flow_all_str[initial_medoids_vec[j]])));
                    if (cal_distance < distance) {
                        distance = cal_distance;
                        cloest_medoids = initial_medoids_vec[j];
                    }      
            }
            closest_medoids_vec.push_back(cloest_medoids);
            absolute_error += distance;
        }
    }
    return absolute_error;
}

// Functions to calculate absolute error for pair of data points by Euclidean Distance
double get_AbsoluteError_E(vector<string> flow_all_str, vector<int> initial_medoids_vec){
    double absolute_error=0;
    vector<int>closest_medoids_vec;
    for(int i=0; i<flow_all_str.size();i++) {
        // excldue initial medoids itselves
        if(contain(initial_medoids_vec, i)){
            closest_medoids_vec.push_back(i);  
        }
            else {
                double distance = 1000000;     // Set an extreme high default value for distance
                int cloest_medoids;
                for(int j=0; j<initial_medoids_vec.size(); j++) {
                        double cal_distance = sqrt(pow(((avg_time(flow_all_str[i])-avg_time(flow_all_str[initial_medoids_vec[j]]))),2.0)
                                            + pow(((avg_length(flow_all_str[i])-avg_length(flow_all_str[initial_medoids_vec[j]]))),2.0));
                        if (cal_distance < distance) {
                            distance = cal_distance;
                            cloest_medoids = initial_medoids_vec[j];
                        }      
                }
                closest_medoids_vec.push_back(cloest_medoids);
                absolute_error += distance;
            }
        }
    return absolute_error;
}

// Functions to get closest mediods list (Manhattan Distance)
vector<int> get_CloestMediodsVec (vector<string> flow_all_str, vector<int>initial_medoids_vec){
    vector<int>closest_medoids_vec;
    double absolute_error =0;
    for(int i=0; i<flow_all_str.size();i++) {
        // excldue initial medoids itselves
        if(contain(initial_medoids_vec, i)){
            closest_medoids_vec.push_back(i);  
        }
            else {
                double distance = 1000000;      // Set an extreme high default value for distance
                int cloest_medoids;
                for(int j=0; j<initial_medoids_vec.size(); j++) {
                        double cal_distance = abs((avg_time(flow_all_str[i])-avg_time(flow_all_str[initial_medoids_vec[j]]))) 
                                            + abs((avg_length(flow_all_str[i])-avg_length(flow_all_str[initial_medoids_vec[j]])));
                        if (cal_distance < distance) {
                            distance = cal_distance;
                            cloest_medoids = initial_medoids_vec[j];
                        }      
                }
                closest_medoids_vec.push_back(cloest_medoids);
                absolute_error += distance;
            }
    }
    return closest_medoids_vec;
}

// Functions to get closest mediods list (Euclidean Distance)
vector<int> get_CloestMediodsVec_E (vector<string> flow_all_str, vector<int>initial_medoids_vec){
    vector<int>closest_medoids_vec;
    double absolute_error =0;
    for(int i=0; i<flow_all_str.size();i++) {
        // excldue initial medoids itselves
        if(contain(initial_medoids_vec, i)){
            closest_medoids_vec.push_back(i);  
        }
            else {
                double distance = 1000000;      // Set an extreme high default value for distance
                int cloest_medoids;
                for(int j=0; j<initial_medoids_vec.size(); j++) {
                        double cal_distance = sqrt(pow(((avg_time(flow_all_str[i])-avg_time(flow_all_str[initial_medoids_vec[j]]))),2.0)
                                            + pow(((avg_length(flow_all_str[i])-avg_length(flow_all_str[initial_medoids_vec[j]]))),2.0));
                        if (cal_distance < distance) {
                            distance = cal_distance;
                            cloest_medoids = initial_medoids_vec[j];
                        }      
                }
                closest_medoids_vec.push_back(cloest_medoids);
                absolute_error += distance;
            }
    }
    return closest_medoids_vec;
}

// Get Distribution 
vector<string> get_Distribution(vector<int>current_medoids_vec, vector<int>final_cloest_mediods){
    vector<string> dist_all;
    for (int i=0; i<current_medoids_vec.size();i++){
        string dist_row="";
        for(int j=0; j<final_cloest_mediods.size();j++) {
            if (current_medoids_vec[i] == final_cloest_mediods[j]) {
                dist_row += to_string(j);
                dist_row += " ";
            }
        }
        dist_all.push_back(dist_row);
    }
    return dist_all;
}

// Get total absolute error (Manhattan Distance)
double get_FinalAbsoluteError(vector<string> flow_all_str, vector<int>initial_medoids_vec){
    // Get the absolute error value for initital mediods
    double current_absolute_error =  get_AbsoluteError(flow_all_str, initial_medoids_vec);
    // Create variable for new abosulte error
    double new_absolute_error=0;
    // Create 2 vectors to store initial mediods
    vector<int>tmp_medoids_vec = initial_medoids_vec;
    vector<int>current_medoids_vec = initial_medoids_vec;

    while(new_absolute_error < current_absolute_error) {
        label1:
        for (int h=0; h<flow_all_str.size(); h++) {
            if((contain(current_medoids_vec, h))==0) {
                for (int i =0; i <current_medoids_vec.size(); i++) {
                    tmp_medoids_vec[i] = h;
                    new_absolute_error = get_AbsoluteError(flow_all_str,tmp_medoids_vec);

                    if(new_absolute_error < current_absolute_error) {
                        // Swapping
                        current_medoids_vec[i] = tmp_medoids_vec[i];    //Update the current mediod
                        current_absolute_error = new_absolute_error;    //Update the current error
                        goto label1;
                    }
                        else {
                            // If haven't do the swapping, change back to the original number
                            tmp_medoids_vec[i] = current_medoids_vec[i];
                        }
                }
            }
        }
    }
    return current_absolute_error;
}

// Get total absolute error (Manhattan Distance)
double get_FinalAbsoluteError_E(vector<string> flow_all_str, vector<int>initial_medoids_vec){
    // Get the absolute error value for initital mediods
    double current_absolute_error =  get_AbsoluteError_E(flow_all_str, initial_medoids_vec);
    // Create variable for new abosulte error
    double new_absolute_error=0;
    // Create 2 vectors to store initial mediods
    vector<int>tmp_medoids_vec = initial_medoids_vec;
    vector<int>current_medoids_vec = initial_medoids_vec;

    while(new_absolute_error < current_absolute_error) {
        label2:
        for (int h=0; h<flow_all_str.size(); h++) {
            if((contain(current_medoids_vec, h))==0) {
                for (int i =0; i <current_medoids_vec.size(); i++) {
                    tmp_medoids_vec[i] = h;
                    new_absolute_error = get_AbsoluteError_E(flow_all_str,tmp_medoids_vec);

                    if(new_absolute_error < current_absolute_error) {
                        // Swapping
                        current_medoids_vec[i] = tmp_medoids_vec[i];    //Update the current mediod
                        current_absolute_error = new_absolute_error;    //Update the current error
                        goto label2;
                    }
                        else {
                            // If haven't do the swapping, change back to the original number
                            tmp_medoids_vec[i] = current_medoids_vec[i];
                        }
                }
            }
        }
    }
    return current_absolute_error;
}

// Get final mediods list (Euclidean Distance)
vector<int> get_FinalMediodsVec(vector<string> flow_all_str, vector<int>initial_medoids_vec){
    // Get the absolute error value for initital mediods
    double current_absolute_error =  get_AbsoluteError(flow_all_str, initial_medoids_vec);
    // Create variable for new abosulte error
    double new_absolute_error=0;
    // Create 2 vectors to store initial mediods
    vector<int>tmp_medoids_vec = initial_medoids_vec;
    vector<int>current_medoids_vec = initial_medoids_vec;

    while(new_absolute_error < current_absolute_error) {
        label3:
        for (int h=0; h<flow_all_str.size(); h++) {
            if((contain(current_medoids_vec, h))==0) {
                for (int i =0; i <current_medoids_vec.size(); i++) {
                    tmp_medoids_vec[i] = h;
                    new_absolute_error = get_AbsoluteError(flow_all_str,tmp_medoids_vec);

                    if(new_absolute_error < current_absolute_error) {
                        // Swapping
                        current_medoids_vec[i] = tmp_medoids_vec[i];    //Update the current mediod
                        current_absolute_error = new_absolute_error;    //Update the current error
                        goto label3;
                    }
                        else {
                            // If haven't do the swapping, change back to the original number
                            tmp_medoids_vec[i] = current_medoids_vec[i];
                        }
                }
            }
        }
    }
    return current_medoids_vec;
}

// Get final mediods list (Euclidean Distance)
vector<int> get_FinalMediodsVec_E(vector<string> flow_all_str, vector<int>initial_medoids_vec){
    // Get the absolute error value for initital mediods
    double current_absolute_error =  get_AbsoluteError_E(flow_all_str, initial_medoids_vec);
    // Create variable for new abosulte error
    double new_absolute_error=0;
    // Create 2 vectors to store initial mediods
    vector<int>tmp_medoids_vec = initial_medoids_vec;
    vector<int>current_medoids_vec = initial_medoids_vec;

    while(new_absolute_error < current_absolute_error) {
        label4:
        for (int h=0; h<flow_all_str.size(); h++) {
            if((contain(current_medoids_vec, h))==0) {
                for (int i =0; i <current_medoids_vec.size(); i++) {
                    tmp_medoids_vec[i] = h;
                    new_absolute_error = get_AbsoluteError_E(flow_all_str,tmp_medoids_vec);

                    if(new_absolute_error < current_absolute_error) {
                        // Swapping
                        current_medoids_vec[i] = tmp_medoids_vec[i];    //Update the current mediod
                        current_absolute_error = new_absolute_error;    //Update the current error
                        goto label4;
                    }
                        else {
                            // If haven't do the swapping, change back to the original number
                            tmp_medoids_vec[i] = current_medoids_vec[i];
                        }
                }
            }
        }
    }
    return current_medoids_vec;
}

// Fuctions to create KMedoids.txt - calculating Error by Manhattan Distance (Question 2)
void Manhattan(vector<string> flow_all_str, vector<int>initial_medoids_vec){
    ofstream newFile2;
    newFile2.open("KMedoids.txt");
    // Absolute Error
    //newFile2 << current_absolute_error << endl;
    newFile2 << get_FinalAbsoluteError(flow_all_str, initial_medoids_vec) << endl;
    // Final Mediods Number
    vector<int> final_mediods_vec = get_FinalMediodsVec(flow_all_str,initial_medoids_vec);
    for (int i=0; i< final_mediods_vec.size();i++) {
        newFile2 << final_mediods_vec[i] << " ";
    }
    // Get final cloest mediods
    vector<int>final_cloest_mediods = get_CloestMediodsVec (flow_all_str, final_mediods_vec);
    newFile2 << endl;
    // Print Final Distribution
    vector<string> dist = get_Distribution(final_mediods_vec, final_cloest_mediods);
    for (int i=0; i<dist.size();i++){
        newFile2 << dist[i] << endl;
    }
    newFile2.close();
}

// Fuctions to create KMedoidsE.txt - calculating Error by Euclidean Distance (Question 3)
void Euclidean(vector<string> flow_all_str, vector<int>initial_medoids_vec){
    ofstream newFile3;
    newFile3.open("KMedoidsE.txt");
    // Absolute Error
    newFile3 << get_FinalAbsoluteError_E(flow_all_str, initial_medoids_vec)  << endl;
    // Final Mediods Number
    vector<int> final_mediods_vec_E = get_FinalMediodsVec_E(flow_all_str, initial_medoids_vec);
    for (int i=0; i< final_mediods_vec_E.size();i++) {
        newFile3 << final_mediods_vec_E[i] << " ";
    }
    // Get final closest mediods
    vector<int> final_cloest_mediods = get_CloestMediodsVec_E(flow_all_str, final_mediods_vec_E);
    newFile3 << endl;
    // Print Final Distribution
    vector<string> dist = get_Distribution(final_mediods_vec_E, final_cloest_mediods);
    for (int i=0; i<dist.size();i++){
        newFile3 << dist[i] << endl;
    }
    newFile3.close();
}



int main() {

	// 1st file
	ifstream file_("file1.txt");

    // skip 1st line
	getline(file_, label_name);

    if(file_.is_open()){
        while(file_ >> src_addr >> src_port >> dst_addr >> dst_port >> protocol >> arrival_time >> packet_length) {
            // Store file data into corresponding vectors (7 columns in total)
            src_addr_vec.push_back(src_addr);
            src_port_vec.push_back(src_port);
            dst_addr_vec.push_back(dst_addr);
            dst_port_vec.push_back(dst_port);
            protocol_vec.push_back(protocol);
            arrival_time_vec.push_back(arrival_time);
            packet_length_vec.push_back(packet_length);
        }
    }

	file_.close();

    int traffic_size = src_addr_vec.size(); // Number of rows for whole data
    vector<int> number_added;               // Store grouped clusters in terms of row number
    vector<string> flow_all_str;            // Store all clusters in group in the form of string

    // Data Structure for storing flows
    for(int i = 0; i < traffic_size; i++) {

        bool flag = false;
        string row_num="";
        
        // check if it is already recorded
        if(std::find(number_added.begin(), number_added.end(), i)!= number_added.end()){
            // if so, then skip
            continue;
        }
            else {
                for (int j=i+1; j<traffic_size; j++) {
                    // Compare values in first 5 columns for all rows
                    if( src_addr_vec[i]==src_addr_vec[j] &&
                        src_port_vec[i]==src_port_vec[j] &&
                        dst_addr_vec[i]==dst_addr_vec[j] &&
                        dst_port_vec[i]==dst_port_vec[j] &&
                        protocol_vec[i]==protocol_vec[j]) {
                        
                        flag = true;    // if flag equals true, we will include the target row
                        row_num += to_string(j);
                        row_num += " ";
                        number_added.push_back(j);
                    }
                }
            }
        if(flag == true) {
            // Add back the target row
            number_added.push_back(i);
            row_num += to_string(i);
            flow_all_str.push_back(row_num);
        } 
    }
    

    // Create Flow.txt
    ofstream newFile;
    newFile.open("Flow.txt");

    for(int i=0; i<flow_all_str.size();i++){
        newFile << i << " ";
        newFile << fixed << setprecision(2) << avg_time(flow_all_str[i]) << " ";
        newFile << fixed << setprecision(2) << avg_length(flow_all_str[i]) << endl;
    }

    newFile.close();

    // 2nd file
	ifstream file2("file2.txt");

    // create variable for storing 1st and 2nd line
    string k_num;
    string second_line;

    // read the 1st line and 2nd line
    getline(file2, k_num);
	getline(file2, second_line);
    second_line = second_line.substr (0,second_line.size()-1);

    // Convert 2nd line(string) into vector of initial medoids(int)
    vector<int>initial_medoids_vec = cluster_int(second_line);

    
    // Question 2 - Create KMedoids.txt
    Manhattan(flow_all_str, initial_medoids_vec);

    // Question 3 - Create KMedoidsE.txt
    Euclidean(flow_all_str, initial_medoids_vec);

    return 0;
}