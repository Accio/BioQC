 #ifndef BIOQC_READ_GMT
 #define BIOQC_READ_GMT

 #include <string>
 #include <vector>
 #include <iostream>
 #include <iterator>
 #include <algorithm>

 class GmtItem {
   private:
   std::string gs_name;
   std::string gs_desc;
   std::vector<std::string> gs_genes;

  public:
   GmtItem(std::istream&);
   GmtItem(const std::string& line);
   GmtItem(std::string const name,
	   std::string const desc,
	   const std::vector<std::string>& genes);

   std::istream& read(std::istream& in);
   void read(const std::string& str);

   std::string name() const {return gs_name;}
   std::string desc() const {return gs_desc;}
   const std::vector<std::string>& genes() const {return gs_genes;} // returns a reference to the vector
   const std::size_t size() const {return gs_genes.size();};
   bool isValid() const {return gs_genes.size()>0;};
 };

 GmtItem::GmtItem(std::istream& in) {read(in);}
 GmtItem::GmtItem(const std::string& line) {read(line);}
 GmtItem::GmtItem(std::string const name, std::string const desc,
		  const std::vector<std::string>& genes) {
   gs_name=name;
   gs_desc=desc;
   copy(genes.begin(), genes.end(), back_inserter(gs_genes));
 }

 bool istab(char c) { return c=='\t'; }
 bool not_tab(char c) { return !istab(c); }

 template <class Out> void split(const std::string& str, Out os) {
   typedef std::string::const_iterator iter;

   iter i=str.begin();
   while(i != str.end()) {
     // ignoring leading blanks
     i = find_if(i, str.end(), not_tab);
     // find end of next word
     iter j=find_if(i, str.end(), istab);

     // copy the characters in [i,j)
     if(i!=str.end())
       *os++=std::string(i,j);

     i=j;
   }
 }

void GmtItem::read(const std::string& s) {
  std::vector<std::string> vec;
  split(s, back_inserter(vec));
  if(vec.size()>=3) { // at least one gene
    std::vector<std::string>::iterator it=vec.begin();
    gs_name=*it++;
    gs_desc=*it++;
    gs_genes.clear();
    copy(it, vec.end(), back_inserter(gs_genes));
  }
}

std::istream& GmtItem::read(std::istream& in) {
  if(in) {
    std::string str;
    if(std::getline(in, str)) {
      read(str);
    }
  }
  return in;
}



#endif
