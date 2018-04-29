#include <fstream>
#include <sstream>
#include <string>
#include <GridPoint.h>

void density_to_vol(grid & points) {
	int NX = points.size();
	int NY = points[0].size();
	int NZ = points[0][0].size();
	int num_channels = 1;
	int xmin = 0
	int ymin = 0
	int zmin = 0
	int xmax = LX * NX;
	int ymax = LY * NY;
	int zmax = LZ * NZ;


	stringstream ss;
	ss<<"VOL";
	ss<<"3";
	ss<<1;
	ss<<NX<<NY<<NZ;
	ss<<num_channels;
	ss<<xmin<<ymin<<zmin<<xmax<<ymax<<zmax;

	for(int i =0 ;i<NX;i++)
	{
		for(int j = 0;j < NY; j++)
		{
			for(int k = 0; k < NZ; k++)
			{
				ss<<points[i][j][k].get_density();
			}
		}
	}
}