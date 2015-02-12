#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>


int main(int argc, char **argv)
{
    FILE *fp = fopen(argv[1], "r");
    assert(fp != NULL);

    // get number of particles
    int np = 0;
    char line[1024];
    while (fgets(line, 1024, fp) != NULL) {
        int nsc = sscanf(line, "%d\n", &np);
        if (nsc != EOF) {
            break;
        }
    }
    assert(np > 0);
    int n3 = np * 3;
    
    // get number of steps
    int nlines = 0;
    rewind(fp);
    double tslen;
    while (fgets(line, 1024, fp) != NULL) {
        double ts0;
        double ts1;
        char strbuf[16];
        int nsc = sscanf(line, "%s %lf\n", strbuf, &ts1);
        if (nsc == EOF) {
            continue;
        }
        if (nlines == 1) {
            ts0 = ts1;
        } else if (nlines == 3 + np) {
            tslen = ts1 - ts0;
        }
        nlines++;
    }
    assert(nlines%(2 + np) == 0);
    int num_steps = nlines/(2 + np);
    printf("nsteps = %d, step-length = %.3le\n", num_steps, tslen);
    assert(num_steps > 2);
    
    // get pos
    std::vector<std::vector<double> > pos (num_steps);
    for (int i = 0; i < num_steps; i++)
    {
        pos[i].reserve(n3);
    }
    int num_radii = 0;
    std::vector<double> radii (np);
    std::vector<int> radius_idx (np);
    std::vector<int> radius_count (np, 0);
    rewind(fp);
    nlines = 0;
    while (fgets(line, 1024, fp) != NULL) {
        char strbuf[16];
        int nsc = sscanf(line, "%s\n", strbuf);
        if (nsc == EOF) {
            continue;
        }
        int step_id = nlines/(2 + np);
        int pos_id = nlines%(2 + np) - 2;
        if (pos_id >= 0) {
            double x;
            double y;
            double z;
            double a;
            nsc = sscanf(line, "%s %lf %lf %lf\n", strbuf, &x, &y, &z);
            a = 1.0;
            pos[step_id][3 * pos_id + 0] = x;
            pos[step_id][3 * pos_id + 1] = y;
            pos[step_id][3 * pos_id + 2] = z;
            if (step_id == 0) {
                int k;
                for (k = 0; k < num_radii; k++) {
                    if (a == radii[k]) {
                        break;    
                    }
                }
                if (k == num_radii) {
                    radii[num_radii] = a;
                    num_radii++;
                }
                radius_idx[pos_id] = k;
                radius_count[k]++;
            }
        }
        nlines++;        
    }
    fclose(fp);
    printf("Have %d different radii\n", num_radii);

    // compute diffusion
    std::vector<double> diffusion (num_radii);
    printf("\ntau, ");
    for (int j = 0; j < num_radii - 1; j++)
    {
        printf("r = %.3f, ", radii[j]);
    }
    printf("r = %.3f\n", radii[num_radii - 1]);
    
    for (int i = 1; i < num_steps/2; i++)
    {
        int count = 0;
        memset(&diffusion[0], 0, num_radii * sizeof (double));
        for (int j = 0; j < num_steps/2; j++)
        {
            int next = j + i;
            for (int k = 0; k < np; k++)
            {
                int idx = radius_idx[k];
                diffusion[idx] +=
                    (pos[next][3 * k] - pos[j][3 * k]) *
                    (pos[next][3 * k] - pos[j][3 * k]);
                diffusion[idx] +=
                    (pos[next][3 * k + 1] - pos[j][3 * k + 1]) *
                    (pos[next][3 * k + 1] - pos[j][3 * k + 1]);
                diffusion[idx] +=
                    (pos[next][3 * k + 2] - pos[j][3 * k + 2]) *
                    (pos[next][3 * k + 2] - pos[j][3 * k + 2]);
            }
            count++;
        }
        for (int j = 0; j < num_radii; j++)
        {
            diffusion[j] /= (6.0 * radius_count[j] * count * i * tslen);
            diffusion[j] *= radii[j];
        }
        
        printf("%.3f, ", i * tslen);
        for (int j = 0; j < num_radii - 1; j++)
        {
            printf("%.6f, %.6f, ", diffusion[j], diffusion[j] * i * tslen);
        }
        printf("%.6f, %.6f\n",
            diffusion[num_radii - 1], diffusion[num_radii - 1] * i * tslen);
    }

        
    return 0;   
}
