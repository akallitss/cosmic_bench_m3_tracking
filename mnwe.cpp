#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

class Base{
	public:
		Base();
		~Base();
		virtual void test() const = 0;
};
class Derived: public Base{
	public:
		Derived();
		~Derived();
		void test() const;
};

Base::Base(){
	cout << "Base constructor called" << endl;
}

Derived::Derived(){
	cout << "Derived constructor called" << endl;
}

Base::~Base(){
	cout << "Base destructor called" << endl;
}

Derived::~Derived(){
	cout << "Derived destructor called" << endl;
}

void Derived::test() const{
	cout << "blah" << endl;
}

int main(){
	vector<Derived> a;
	cout << 1 << endl;
	Derived b;
	cout << 2 << endl;
	a.push_back(b);
	cout << 3 << endl;
	vector<Derived>::iterator jt = a.end();
	for(vector<Derived>::iterator it = a.begin();it!=a.end();++it){
		jt = it;
	}
	jt->test();
	return 0;
}