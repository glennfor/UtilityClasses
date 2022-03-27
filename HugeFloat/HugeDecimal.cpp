using namespace std;

#if 0
//all operators
//arithmetic
 + - / * % and - to negate

//bitwise
^ & | >>  <<  ~//might use shift operators for stream insertions, no need

//boolean
|| && ! > < == >=  <=

//increment and decrement
++ -- //(same implemntaion for pre and post)


//
 
 #endif
class BigInteger{
	friend operator+(BigInteger&, BigInteger&);
	friend operator+(BigInteger&, long);//will also be able to take int
	
	friend operator*(BigInteger&, BigInteger&);
	friend operator*(BigInteger&, long);
	
	friend operator/(BigInteger&, BigInteger&);
	friend operator/(BigInteger&, long);
	
	friend operator-(BigInteger&, BigInteger&);
	friend operator-(BigInteger&, long);
	
	private:
		std:: string number_data;
	public:
		BigInteger(std::string num):number_data(num){}
		BigInteger(long long number){
			//use stringstream or sprintf to insert number
		}
		std::string toString()
		{
			return number_data;
		}
		BigInteger operator=(BigInteger& bigint)
		{
			this->number_data = bigint.toString();
			return *this;
		}
 		BigInteger operator=(long long longint)
		{
			//use stringstream
			this->number_data = 
			return *this;
		}
};

