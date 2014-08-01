/**
 * @file   molecule_io.cc
 * @brief  the MoleculeIO class implementation
 */

#include <cstdio>
#include <cstring>
#include <iterator>
#include "molecule_io.h"
#include "log.h"


namespace stokesdt{

using namespace std;
    
/** @brief  Construct a MoleculeIO object
 */
MoleculeIO::MoleculeIO()
{

}


/** @brief  Deconstruct MoleculeIO
 */
MoleculeIO::~MoleculeIO()
{

}


/** @brief  Import keys and values from the config file
 *
 *  @param[in] config_file  String containing the config file name
 *  @param[in] deleimters   String containing the delimiter characters
 *  @param[in] comments     String containing the prefix of comment lines
 */
bool MoleculeIO::ParseConfig(const char *config_file,
                             const char *delimiters,
                             const char *comments)
{
    FILE *fp = fopen(config_file, "r");
    if (NULL == fp) {
        LOG_ERROR("Open file %s failed\n", config_file);
    }
    char line[detail::kMaxLine];
    string str_key;
    string str_value;
    key_map_.clear();
    while (fgets(line, detail::kMaxLine, fp) != NULL) {
        char *line0 = strdup(line);
        char *head = strtok(line0, " \t");
        // not a comment line
        if (strlen(head) < strlen(comments) ||
            strncmp(head, comments, strlen(comments)) != 0) {
            char *key = strtok(line, delimiters);
            char *value = strtok(NULL, delimiters);
            // a valid line
            if (key != NULL && value != NULL) {
                str_key.assign(key);
                str_value.assign(value);
                if (key_map_.find(str_key) == key_map_.end()) {
                    key_map_[str_key] = str_value;    
                } else {
                    LOG_ERROR("Duplicate keys (%s) found in %s\n",
                        key, config_file);
                    return false;
                }
            }
        }
        free(line0);
    }
    
    fclose(fp);
    return true;
}


/** @brief  Import model lines from the model file
 *
 *  @param[in] config_file  String containing the model file name
 *  @param[in] deleimters   String containing the delimiter characters
 *  @param[in] comments     String containing the prefix of comment lines
 */
bool MoleculeIO::ParseModel(const char *model_file,
                            const char *delimiters,
                            const char *comments)
{
    FILE *fp = fopen(model_file, "r");
    if (NULL == fp) {
        LOG_ERROR("Open file %s failed\n", model_file);
    }
    char line[detail::kMaxLine];
    string str_key;
    string str_value;
    line_map_.clear();
    while (fgets(line, detail::kMaxLine, fp) != NULL) {
        char *line0 = strdup(line);
        char *head = strtok(line0, " \t");
        // not a comment line
        if (strlen(head) < strlen(comments) ||
            strncmp(head, comments, strlen(comments)) != 0) {
            char *key = strtok(line, delimiters);
            char *value = strtok(NULL, "");
            // a valid line
            if (key != NULL && value != NULL) {
                str_key.assign(key);
                str_value.assign(value);
                line_map_.insert(pair<string, string>(str_key, str_value));
            }
        }
        free(line0);
    }
    
    fclose(fp);
    return true;
}


/** @brief  Import particle radii and positions
 *
 *  @param[in] xyz_file     String containing the XYZ file name
 *  @param[out] npos        Return the number of particles
 *  @param[out] Lx          Return dimension along x-axis
 *  @param[out] Ly          Return dimension along y-axis
 *  @param[out] Lz          Return dimension along z-axis
 *  @param[out] start_step  Return step index
 *  @param[out] start_time  Return start time
 */
bool MoleculeIO::ParseXYZ(const char *xyz_file,
                          int *npos, double *Lx,
                          double *Ly, double *Lz,
                          int *start_step, double *start_time)
{
    // read particle types
    int num_types = LineCount("par");
    if (num_types == 0) {
        LOG_WARN("Paticle types not found. Set radii to 1.0\n");
    }
    map<string, double> par_map;
    string str_line;
    string str_par;
    char type[128];
    for (int i = 0; i < num_types; i++) {
        double radius;
        str_line = GetLine("par", i);
        sscanf(str_line.c_str(), "%s %le", type, &radius);
        str_par.assign(type);
        par_map[str_par] = radius;
    }

    // read XYZ
    FILE *fp = fopen(xyz_file, "r");
    if (NULL == fp) {
        LOG_ERROR("Open file %s failed\n", xyz_file);
    }
    char line[detail::kMaxLine];
    // read npos
    if (fgets(line, detail::kMaxLine, fp) != NULL) {
        *npos = atoi(line);
        if (*npos <= 0) {
            LOG_ERROR("Invalid XYZ file %s\n", xyz_file);
            return false;    
        }
    } else {
        LOG_ERROR("Invalid XYZ file %s\n", xyz_file);
        return false;
    }
    // read the comment line
    if (fgets(line, detail::kMaxLine, fp) != NULL) {
        sscanf(line, "%d %le %le %le %le",
            start_step, start_time, Lx, Ly, Lz);
        if (*start_step < 0 || *start_time < 0.0 ||
            *Lx <= 0.0 || *Ly <= 0.0 || *Lz <= 0.0) {
            LOG_ERROR("Invalid XYZ file %s\n", xyz_file);
            return false;            
        }
    } else {
        LOG_ERROR("Invalid XYZ file %s\n", xyz_file);
        return false;
    }    
    // read pos and rdi 
    pos_.clear();
    rdi_.clear(); 
    while (fgets(line, detail::kMaxLine, fp) != NULL) {
        double x;
        double y;
        double z;
        sscanf(line, "%s %le %le %le", type, &x, &y, &z);
        if (strtok(line, " \t") == NULL) {
            continue;
        }
        str_par.assign(type);
        pos_.push_back(x);
        pos_.push_back(y);
        pos_.push_back(z);
        if (par_map.find(str_par) != par_map.end()) {
            rdi_.push_back(par_map[str_par]);
        } else {
            rdi_.push_back(1.0);
        }   
    }
    
    // check npos
    if (*npos != rdi_.size()) {
        LOG_ERROR("Invalid XYZ file %s\n", xyz_file);
        return false;
    }

    fclose(fp);
    return true;
}


/** @brief  Import particle radii and positions
 *
 *  @param[in] xyz_file     String containing the XYZ file name
 *  @param[out] npos        Return the number of particles
 *  @param[out] box_size    Return dimension of simulation box
 */
bool MoleculeIO::ParseXYZ(const char *xyz_file, int *npos, double *box_size)
{
    double Ly;
    double Lz;
    int start_step;
    double start_time;
    
    if (ParseXYZ(xyz_file, npos, box_size, 
            &Ly, &Lz, &start_step, &start_time)) {
        return true;
    } else {
        return false;
    }
}


double MoleculeIO::GetFloatKey(const char *key,
                               const double min_value,
                               const double max_value,
                               const double default_value)
{
    string str_key(key);
    
    if (key_map_.find(str_key) != key_map_.end()) { // key not found
        LOG_WARN("Key \"%s\" not found. Use default value: %d\n",
                 key, default_value);
        return default_value;
    } else {
        double ret = atof(key_map_[str_key].c_str());
        if (ret >= max_value || ret <= min_value) {
            LOG_WARN("Key \"%s\" (%lf) not in the range (%lf %lf). "
                     "Use default value: %lf\n",
                     key, ret, min_value, max_value, default_value);
            return default_value;
        } else {
            return ret;
        }
    }
}


int MoleculeIO::GetIntKey(const char *key,
                          const int min_value,
                          const int max_value,
                          const int default_value)
{
    string str_key(key);
    
    if (key_map_.find(str_key) != key_map_.end()) { // key not found
        LOG_WARN("Key \"%s\" not found. Use default value: %d\n",
                 key, default_value);
        return default_value;
    } else {
        int ret = atoi(key_map_[str_key].c_str());
        if (ret > max_value || ret < min_value) {
            LOG_WARN("Key \"%s\" (%lf) not in the range [%d %d]. "
                     "Use default value: %d\n",
                     key, ret, min_value, max_value, default_value);
            return default_value;
        } else {
            return ret;
        }
    }
}


bool MoleculeIO::GetBoolKey(const char *key,
                            const bool default_value)
{
    string str_key(key);

    if (key_map_.find(str_key) != key_map_.end()) { // key not found
        LOG_WARN("Key \"%s\" not found. Use default value: %s\n",
                 key, (default_value ? "true" : "false"));
        return default_value;        
        return false;
    } else {
        if (key_map_[str_key].compare("true") == 0 ||
            key_map_[str_key].compare("True") == 0 ||
            key_map_[str_key].compare("1") == 0) {
            return true;
        } else if (key_map_[str_key].compare("false") == 0 ||
                   key_map_[str_key].compare("False") == 0 ||
                   key_map_[str_key].compare("0") == 0) {
            return false;
        } else {
             LOG_WARN("Key \"%s\" (%s) not a bool value. "
                     "Use default value: %s\n",
                     key, str_key.c_str(), (default_value ? "true" : "false"));
            return default_value;       
        }
    }
}


const char* MoleculeIO::GetStringKey(const char *key)
{
    string str_key(key);
    
    if (key_map_.find(str_key) != key_map_.end()) { // key not found
        return NULL;
    } else {
        return key_map_[str_key].c_str();
    }
}


int MoleculeIO::LineCount(const char *key)
{
    string str_key(key);
    return line_map_.count(str_key);    
}


const char* MoleculeIO::GetLine(const char *key, int idx)
{
    auto ret = line_map_.equal_range(key);
    multimap<string, string>::iterator it;
    it = ret.first;
    advance(it, idx);
    return it->second.c_str();
}


void MoleculeIO::GetParticles(double *pos, double *rdi)
{
    int npos = rdi_.size();
    memcpy(pos, &pos_[0], sizeof(double) * npos * 3);
    memcpy(rdi, &rdi_[0], sizeof(double) * npos);
}

} // namespace stokesdt
