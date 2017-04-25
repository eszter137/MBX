#include "read_nrg.h" 

namespace tools {

void ReadNrg(char * filename, std::vector<bblock::System> & systems ) {
  assert(filename);
  std::ifstream ifs(filename);

  if (!ifs)
    throw std::runtime_error("could not open NRG file for reading");
  
  size_t lineno(0);
  size_t sysno(0);

  while (true) {
    bblock::System sys;
    ReadSystem(lineno, ifs, sys);
    
    systems.push_back(sys);
    sysno++;

    std::streampos oldpos = ifs.tellg();
    std::string line;
    std::getline(ifs, line);
    if (ifs.eof() or line.empty()) {
      break;
    } else {
      ifs.seekg(oldpos);
    }
  }
  return;
}

void ReadSystem(size_t& lineno, std::istream& ifs, bblock::System& sys) {
  assert(ifs);

  if (ifs.eof())
    return;

  std::string line;
  std::getline(ifs, line);
  lineno++;

  std::string word;
  std::istringstream iss(line);
  iss >> word;
  
  if (iss.fail()) {
    std::ostringstream oss;
    oss << "unexpected text at line " << lineno
        << " of the NRG file:" << std::endl << iss.str() << std::endl;
    throw std::runtime_error(oss.str());
  }

  std::transform(word.begin(), word.end(), word.begin(), ::tolower);

  if (word != "system") {
    std::ostringstream oss;
    oss << "No SYSTEM in line  " << lineno
        << " of the NRG file:" << std::endl << iss.str() << std::endl; 
    throw std::runtime_error(oss.str());
  }

  size_t molno(0);

  while (word != "endsys") {
    std::shared_ptr<bblock::Molecule> molec = 
        std::shared_ptr<bblock::Molecule> (new bblock::Molecule);
    ReadMolecule(lineno, ifs, molec);

    if (molec->GetNumMon() > 0) {
      sys.AddMolecule(molec);
      molno++;
    }
    iss.clear();

    std::streampos oldpos = ifs.tellg();

    std::getline(ifs, line);
    lineno++;
    iss.str(line);
    iss >> word;
    std::transform(word.begin(), word.end(), word.begin(), ::tolower);

    if (word == "molecule") {
      ifs.seekg(oldpos);
      lineno--;
    }
  }

  sys.SetNumMol(molno);
  return;
}

void ReadMolecule(size_t& lineno, std::istream& ifs, std::shared_ptr<bblock::Molecule> molec) {
  assert(ifs);

  if (ifs.eof())
    return;

  std::string line;
  std::getline(ifs, line);
  lineno++;

  std::string word;
  std::istringstream iss(line);
  iss >> word;
  
  if (iss.fail()) {
    std::ostringstream oss;
    oss << "unexpected text at line " << lineno
        << " of the NRG file:" << std::endl << iss.str() << std::endl;
    throw std::runtime_error(oss.str());
  }

  std::transform(word.begin(), word.end(), word.begin(), ::tolower);
  if (word != "molecule") {
    std::ostringstream oss;
    oss << "No MOLECULE in line  " << lineno
        << " of the NRG file:" << std::endl;
    throw std::runtime_error(oss.str());
  }
  
  ReadMonomers(lineno, ifs, molec);

  return;

}

void ReadMonomers(size_t& lineno, std::istream& ifs, std::shared_ptr<bblock::Molecule> molec) {
  assert(ifs);

  if (ifs.eof())
    return;

  std::string line;
  std::getline(ifs, line);
  lineno++;

  std::string word;
  std::istringstream iss(line);
  iss >> word;

  if (iss.fail()) {
    std::ostringstream oss;
    oss << "unexpected text at line " << lineno
        << " of the NRG file:" << std::endl << iss.str() << std::endl;
    throw std::runtime_error(oss.str());
  }

  std::transform(word.begin(), word.end(), word.begin(), ::tolower);
  if (word != "monomer") {
    std::ostringstream oss;
    oss << "No MONOMER in line  " << lineno
        << " of the NRG file:" << std::endl;
    throw std::runtime_error(oss.str());
  }

  std::string mon_name;
  iss.clear();
  iss.str(line);
  iss >> word >> mon_name;
  std::transform(word.begin(), word.end(), word.begin(), ::tolower);

  size_t monno(0);

  while (word != "endmol") {
  
    if (word != "monomer" || mon_name.empty() ) {
      std::ostringstream oss;
      oss << "No MONOMER in line  " << lineno
          << " of the NRG file:" << std::endl << iss.str() << std::endl;
      throw std::runtime_error(oss.str());
    }
    
    std::transform(mon_name.begin(), mon_name.end(), mon_name.begin(), ::tolower);

    std::vector<std::string> names;
    double xyz[10000];
    size_t i(0);
    
    std::getline(ifs, line);
    lineno++;
    iss.clear();
    iss.str(line);
 
    while (word != "endmon") {

      std::string name;
      double x,y,z;
      
      iss >> name >> x >> y >> z;
      names.push_back(name);
      xyz[i++] = x;
      xyz[i++] = y;
      xyz[i++] = z;

      std::getline(ifs, line);
      lineno++;
      iss.clear();
      iss.str(line);

      iss >> word;
      std::transform(word.begin(), word.end(), word.begin(), ::tolower);
      iss.clear();
      iss.str(line);
    }

    molec->AddMonomer(mon_name, xyz, names);
    monno++;

    iss.clear();
    std::getline(ifs, line);
    lineno++;
    iss.str(line);
    iss >> word;
    std::transform(word.begin(), word.end(), word.begin(), ::tolower);
  }
  
  molec->SetNumMon(monno);

  return;
}



} // namespace tools
















