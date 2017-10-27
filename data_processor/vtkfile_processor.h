#ifndef __VTKFILE_PROCESSOR__
#define __VTKFILE_PROCESSOR__

#include <string>
#include <vector>
class NeighbourSearcher;

class VTKFileProcessor {
public:
	VTKFileProcessor();
	~VTKFileProcessor();
	void read(const std::string& filename);
	void write(const std::string& filename);
	void rotate(double angle);
	void translate(double xshift, double yshift, double zshift);
	void combine(const std::string& filename);
	void changeTag(int tag);
	void changeVelocity(double u_, double v_, double w_);
	void reduce2dY(double yCen, double xLeft, double xRight, double spacing, const std::string& fname);
	void reduce3dYZ(double yCen, double Zcen, double xLeft, double xRight, double spacing, const std::string& fname);
	static std::string rightFlush(std::size_t writeStep, std::size_t numDigits);
        void interpolation_to_mesh(double xLeft, double xRight, double yLeft, double yRight, double zLeft, double zRight, int xN, int yN, int zN, const std::string& filename);
	void writeMesh(double* massMesh, double xLeft, double xRight, double yLeft, double yRight, double zLeft, double zRight, int xN, int yN, int zN, const std::string& filename);
private:
	std::size_t numParticle;
	double time;
	NeighbourSearcher* ns;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<double> u;
	std::vector<double> v;
	std::vector<double> w;
	std::vector<double> volume;
	std::vector<double> pressure;
	std::vector<double> soundSpeed;
	std::vector<int> objectTag;
	std::vector<double> localParSpacing;
	std::vector<double> mass;

	void alloc();
	void dealloc();
	void augment(std::size_t added_size);
	double linear_interpolate(double diff0, double value0, double diff1, double value1);
	void write1d(const std::string& fname, const std::vector<std::string>& titles, 
	const std::vector<double*>& inputs, std::size_t num);
};

#endif // __VTKFILE_PROCESSOR__
