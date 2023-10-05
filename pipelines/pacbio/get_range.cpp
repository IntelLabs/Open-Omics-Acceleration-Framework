#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

class meta {
	public:
		string cand_file_name;
		string region_file_name;
		int64_t start;
		int64_t end;
		int64_t size;
		meta(string cf, string rf, int64_t sz){
			cand_file_name = cf;
			region_file_name = rf;
			size = sz;
			start = 0;
			end = size - 1;
		}
		meta(){

		}
		int64_t get_size(){
			return end - start + 1;
		}

		void print(){
			cout<<cand_file_name<< " " <<region_file_name << " " << start << " " << end<<endl;
		}
};

vector<meta> output;
vector<vector<meta> > l_output;

void generate_output(meta* m, int64_t end = -1){
	meta ret_m;
	
	ret_m.cand_file_name = m->cand_file_name;
	ret_m.region_file_name = m->region_file_name;
	ret_m.start = m->start;

	ret_m.end = (end == -1) ? m->end: end;

	ret_m.size = ret_m.get_size();

	output.push_back(ret_m);

	//ret_m.print();

}


void fill_buffer(meta* m, int64_t &buff_size, string &tasks) {

	// 1. if size is smaller than buffer_size, fill buffer
	if(m->get_size() <= buff_size) {
		// output
		generate_output(m);
	

		if(m->get_size() == buff_size){
			l_output.push_back(output);
			output.clear();
		}
		buff_size -= m->get_size();
		m->start = m->end + 1;
		return;
	}
	// 2. if size is larger than buffer size
	if ( m->get_size() > buff_size){
		int64_t range_end = ((m->start + buff_size) - 1);
		generate_output(m, range_end);
		m->start = range_end + 1;
		buff_size = 0;
		l_output.push_back(output);
		output.clear();
		return;
	}	
}


string get_range(vector<meta> v, int64_t max_size, int64_t rem){

	string tasks="";
	int i = 0;
	int64_t buff_size;
	if(rem > 0){
		buff_size = max_size + 1;
	
	}
	else {
		buff_size = max_size;
	}
	while(i < v.size()){

	//	v[i].print();
		if(v[i].get_size() > 0){
			fill_buffer(&v[i], buff_size, tasks);
			if(buff_size == 0) {
				buff_size = max_size;
				if(l_output.size() < rem){
					buff_size++;
				}
			}
		}
		else {
			i++;
		}
	}

	if(output.size() != 0){
		l_output.push_back(output);
	}
	return tasks;
}

int main(int argc, char* argv[]) {

	ifstream fp(argv[1]);
	int64_t max_size = atoi(argv[2]);
	
	vector<meta> v;
	
	// File reading
	string cf;
	string rf;
	int64_t cnt;
	int64_t total_count = 0;

	while(fp>>cf){
		fp>>rf;
		fp>>cnt;
		total_count += cnt;
		//cout<<cf<< " " << rf<< " " << cnt <<endl;
		v.push_back(meta(cf, rf, cnt));
	}

	int64_t rem = total_count % v.size();

	fp.close();

	ofstream fp_cand(argv[3]);
	ofstream fp_reg(argv[4]);
	ofstream fp_offset(argv[5]);
	ofstream fp_count(argv[6]);


	cout<<"SMK read records: "<<v.size()<<endl;
	string regions_list = get_range(v, max_size, rem);
	cout<<regions_list;

	for(int i = 0; i < l_output.size(); i++) {

		int count = 0;
		for( int j = 0; j < l_output[i].size(); j++) {
			//l_output[i][j].print();
			fp_cand<<l_output[i][j].cand_file_name << " ";
		    	fp_reg<<l_output[i][j].region_file_name << " ";
			fp_offset<<l_output[i][j].start<<" ";
			count += l_output[i][j].size;	
		}
		fp_count<<count<<endl;
		fp_cand<<"\n";
		fp_reg<<"\n";
		fp_offset<<"\n";	

		//cout<<endl;
	}
		fp_cand.close();
		fp_reg.close();
		fp_offset.close();	
		
}
