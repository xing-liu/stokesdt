/**
 * @file   molecule_io.h
 * @brief  the MoleculeIO class definition
 */
 
#ifndef MOLECULE_IO_H_
#define MOLECULE_IO_H_


#include <cmath>
#include <string>
#include <vector>
#include <map>
#include "common.h"


namespace stokesdt{

/** 
 * @class  MoleculeIO
 * @brief  Wraps the routines for reading and writing molecule files
 */
class MoleculeIO {
  public:
    /// the maximum length the type string
    enum { kMaxTypeLen = 16 };
        
    /// Class constructor
    MoleculeIO();

    /// class deconstructor
    virtual ~MoleculeIO();

    /// Imports config file
    bool ParseConfig(const char *config_file,
                     const char *delimiters,
                     const char *comments);

    /// Imports model file
    bool ParseModel(const char *model_file,
                    const char *delimiters,
                    const char *comments);

    /// Imports XYZ file returning the starting frame and starting time 
    bool ParseXYZ(const char *xyz_file, int *npos, 
                  double *Lx, double *Ly, double *Lz,
                  int *start_frame, double *start_time);

    /// Imports XYZ file
    bool ParseXYZ(const char *xyz_file, int *npos, double *box_size);

    /// Reads a float key
    double GetFloatKey(const char *key, const double min_value,
                       const double max_value, const double default_value);

    /// Reads a integer key
    int GetIntKey(const char *key, const int min_value,
                  const int max_value, const int default_value);

    /// Reads a bool key
    bool GetBoolKey(const char *key, const bool default_value);

    /// Reads a string key
    const char *GetStringKey(const char *key);

    /// Returns the number of specific lines
    int LineCount(const char *key);

    /// Returns a specific line
    const char *GetLine(const char *key, int idx);

    /// Reads particles
    void GetParticles(double *pos, double *rdi);

    /// Reads particles returning ptype
    void GetParticles(double *pos, double *rdi,
                      char *ptype[kMaxTypeLen]);

    /// Writes a XYZ file specifying start_frame, start_time
    void WriteXYZ(const int npos,
                  const double Lx,
                  const double Ly,
                  const double Lz,
                  const int start_frame,
                  const double start_time,
                  const char *ptype[kMaxTypeLen],
                  const double *pos,
                  const double *rdi,
                  FILE* fp);

    /// Export a XYZ file
    void WriteXYZ(const int npos,
                  const double box_size,
                  const char *ptype[kMaxTypeLen],
                  const double *pos,
                  const double *rdi,
                  FILE* fp);
   
  private:
    DISALLOW_COPY_AND_ASSIGN(MoleculeIO);

  private:    
    /// the map storing keys and values
    std::map<std::string, std::string> key_map_;
    /// the map storing lines
    std::multimap<std::string, std::string> line_map_;
    /// the array of particle coordinates
    std::vector<double> pos_;
    /// the array of particle types
    std::vector<std::string> ptype_;
    /// the array of particle radii
    std::vector<double> rdi_;
};

} // namespace stokesdt


#endif // MOLECULE_IO_H_