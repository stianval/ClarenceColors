#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <set>
#include <iterator>
#include <png.h>
using namespace std;

struct Bob {
	string name;
	int index;
	vector<int> colors;
	vector<Bob*> appearsWith;
};

struct Site {
	string name;
	int index;
	vector<int> colors;
	vector<Bob*> bobs;
	vector<size_t> bobEndIndex;

	Bob *currentBob(void) { return bobs.back(); }
};

struct ColorDist {
	size_t indexA, indexB;
	int colorA, colorB;
	int distance;
};

bool operator<(const ColorDist &lhs, const ColorDist &rhs) {
	return lhs.distance < rhs.distance;
}

int hex2dec(const string &hex) {
	int res = 0;
	string::const_iterator it = hex.cbegin() + 1;
	for (; it < hex.cend(); ++it) {
		res *= 16;
		if (*it >= 'a') {
			res += *it - 'a' + 10;
		}
		else if (*it >= 'A') {
			res += *it - 'A' + 10;
		}
		else {
			res += *it - '0';
		}
	}
	return res;
}

void readfile(istream &file, vector<Site> &sites, vector<Bob> &bobs) {
	string line;
	getline(file, line);
	{
		istringstream iss(line);
		string cell;
		while (getline(iss, cell, ','), iss) {
			// Make a new site
			Site site;
			site.name = cell;
			sites.push_back(site);

			// Make a default "background bob" for the site
			Bob bob;
			bob.name = cell + "[bg]";
			bobs.push_back(bob);
			sites.back().bobs.push_back(&bobs.back());
			sites.back().bobEndIndex.push_back(0);
		}
	}
	
	while (getline(file, line), file) {
		istringstream iss(line);
		vector<Site>::iterator siteIt = sites.begin();
	
		for (; siteIt < sites.end(); ++siteIt) {	
			string cell;
			getline(iss, cell, ',');

			// if the cell is empty, skip it.
			if (cell == "") continue;

			// if the cell contains a color (i.e. first letter is #), store it to the site and the appropriate bob.
			if (cell[0] == '#') {
				int color = hex2dec(cell);
				siteIt->colors.push_back(color);
				siteIt->currentBob()->colors.push_back(color);
				siteIt->bobEndIndex.back()++;
			}
			// otherwise the cell contains the name of a bob
			else {
				Bob bob;
				bob.name = cell;
				bobs.push_back(bob);
				siteIt->bobs.push_back(&bobs.back());
				siteIt->bobEndIndex.push_back(siteIt->colors.size());
			}
		}
	}
}

int halveColor(int color) {
	int red = color >> 16;
	int green = color >> 8 & 0xFF;
	int blue = color & 0xFF;
	return red / 2 << 16 | green / 2 << 8 | blue / 2;
}

int amiColor(int color) {
	// extract components
	int red = color >> 16;
	int green = color >> 8 & 0xFF;
	int blue = color & 0xFF;
	
	// convert to amiga color space
	red = red / 0x10;
	green = green / 0x10;
	blue = blue / 0x10;

	// return components combined
	return red * 0x100 | green * 0x10 | blue * 0x1;
}

int squaredEuclideanDistance(int colorA, int colorB) {
	// components of colorA
	int redA = colorA >> 16;
	int greenA = colorA >> 8 & 0xFF;
	int blueA = colorA & 0xFF;

	// components of colorB
	int redB = colorB >> 16;
	int greenB = colorB >> 8 & 0xFF;
	int blueB = colorB & 0xFF;

	// distance from components
	int dRed = redB - redA;
	int dGreen = greenB - greenA;
	int dBlue = blueB - blueA;
	return dRed*dRed + dGreen*dGreen + dBlue*dBlue;
}

int amiSquaredEuclideanDistance(int colorA, int colorB) {
	// components of colorA
	int redA = colorA >> 8;
	int greenA = colorA >> 4 & 0xF;
	int blueA = colorA & 0xF;

	// components of colorB
	int redB = colorB >> 8;
	int greenB = colorB >> 4 & 0xF;
	int blueB = colorB & 0xF;

	// distance from components
	int dRed = redB - redA;
	int dGreen = greenB - greenA;
	int dBlue = blueB - blueA;
	return dRed*dRed + dGreen*dGreen + dBlue*dBlue;
}

void outputDistances(ostream &stream, const vector<int> &allColors, const vector<vector<int> > &distances) {
	stream << setw(6) << "" << " ";
	for (int color : allColors) {
		stream << setw(6) << hex << color << " ";
	}
	stream << '\n';
	for (size_t row = 0; row < distances.size(); ++row) {
		stream << setw(6) << hex << allColors[row] << " ";
		for (size_t col = 0; col <= row; ++col) {
			stream << setw(6) << "" << " ";
		}
		for (int dist : distances[row]) {
			stream << setw(6) << dec << dist << " ";
		}
		stream << '\n';
	}
	stream << flush;
}

bool isHalfbrite(size_t index) {
	return (index & 1) == 1;
}

bool isHalfbritePair(size_t indexA, size_t indexB) {
	return (indexA & 1 ^ indexB & 1) == 1;
}

bool bothHalfbrite(size_t indexA, size_t indexB) {
	return (indexA & indexB & 1) == 1;
}

void outputLeastDistances(ostream &stream, const vector<ColorDist> &sortedDistances, size_t nDistances = 0, int maxAmiDistance = -1, bool onlyHalfbrite = false) {
	if (nDistances == 0) {
		// output all distances
		nDistances = sortedDistances.size();
	}
	if (maxAmiDistance == -1) {
		maxAmiDistance = INT_MAX;
	}
	for (size_t i = 0; i < nDistances; ++i) {
		const ColorDist &colorDist = sortedDistances[i];

		if (bothHalfbrite(colorDist.indexA, colorDist.indexB)) continue;
		if (onlyHalfbrite && !isHalfbritePair(colorDist.indexA, colorDist.indexB)) continue;

		int amiDistance = amiSquaredEuclideanDistance(amiColor(colorDist.colorA), amiColor(colorDist.colorB));

		if (amiDistance > maxAmiDistance) continue;

		// distance data for pc
		bool halfbritePair = isHalfbritePair(colorDist.indexA, colorDist.indexB);
		stream << setw(5) << boolalpha << halfbritePair << " - ";
		stream << '(' << setw(3) << dec << colorDist.indexA << " - " << setw(6) << hex << colorDist.colorA << ')' << ' ';
		stream << '(' << setw(3) << dec << colorDist.indexB << " - " << setw(6) << hex << colorDist.colorB << ')' << ' ';
		stream << setw(6) << dec << colorDist.distance;

		// distance data for amiga
		stream << "    =>    ";
		stream << '(' << setw(3) << dec << colorDist.indexA << " - " << setw(3) << hex << amiColor(colorDist.colorA) << ')' << ' ';
		stream << '(' << setw(3) << dec << colorDist.indexB << " - " << setw(3) << hex << amiColor(colorDist.colorB) << ')' << ' ';
		stream << setw(6) << dec << amiDistance;

		stream << '\n';
	}
	stream << flush;
}

void extractColorComponents(int color, int &red, int &green, int &blue) {
	red = color >> 16;
	green = color >> 8 & 0xFF;
	blue = color & 0xFF;
}

struct Cluster {
	vector<size_t> colorIndices;
	vector<size_t> halfbriteIndices;
	double redMean = 0.0;
	double greenMean = 0.0;
	double blueMean = 0.0;

	void insertColor(int color, size_t index) {
		size_t n = colorIndices.size();
		int red, green, blue;
		extractColorComponents(color, red, green, blue);

		redMean = (redMean * n + red) / (n + 1);
		greenMean = (greenMean * n + green) / (n + 1);
		blueMean = (blueMean*n + blue) / (n + 1);
		colorIndices.push_back(index);
	}

	void insertHalfbrite(int color, size_t index) {
		size_t n = colorIndices.size();
		int red, green, blue;
		extractColorComponents(color, red, green, blue);

		redMean = (redMean * n + 2*red) / (n + 1);
		greenMean = (greenMean * n + 2*green) / (n + 1);
		blueMean = (blueMean*n + 2*blue) / (n + 1);
		halfbriteIndices.push_back(index);
	}

	double distance(int color) {
		int red, green, blue;
		extractColorComponents(color, red, green, blue);

		double dRed = red - redMean;
		double dGreen = green - greenMean;
		double dBlue = blue - blueMean;

		return dRed*dRed + dGreen*dGreen + dBlue*dBlue;
	}

	double halfbriteDistance(int color) {
		int red, green, blue;
		extractColorComponents(color, red, green, blue);

		double dRed = red - redMean / 2.0;
		double dGreen = green - greenMean / 2.0;
		double dBlue = blue - blueMean / 2.0;

		return dRed*dRed + dGreen*dGreen + dBlue*dBlue;
	}
};

void epsilonClusters(vector<Cluster> &result, const vector<int> &allColors, double epsilon) {
	for (size_t index = 0; index < allColors.size(); ++index) {
		bool foundCluster = false;
		for (Cluster &cluster : result) {
			if (cluster.distance(allColors[index]) <= epsilon) {
				cluster.insertColor(allColors[index], index);
				foundCluster = true;
				break;
			}
			if (cluster.halfbriteDistance(allColors[index]) <= epsilon) {
				cluster.insertHalfbrite(allColors[index], index);
				foundCluster = true;
				break;
			}
		}
		if (!foundCluster) {
			Cluster cluster;
			cluster.insertColor(allColors[index], index);
			result.push_back(cluster);
		}
	}
}

void outputBobNames(ostream &stream, const vector<size_t> colorIndices, const vector<size_t> indexMap, const Site &site) {
	set<string> bobNames;
	transform(colorIndices.begin(), colorIndices.end(), inserter(bobNames, bobNames.begin()), [&indexMap, &site](size_t index) {
		size_t mappedIndex = indexMap[index];
		size_t bobNo = 0;
		for (; mappedIndex >= site.bobEndIndex[bobNo]; ++bobNo) { ; }
		string name = site.bobs[bobNo]->name;
		return name;
	});
	for (const string &name : bobNames) {
		stream << name << ' ';
	}
}

void outputClusters(ostream &stream, const vector<Cluster> &clusters, const vector<size_t> &indexMap, const Site &site) {
	int counter = 0;
	stream << "Colors\n";
	for (const Cluster &cluster : clusters) {
		stream << counter++ << " - " << cluster.redMean << ", " << cluster.greenMean << ", " << cluster.blueMean << "   -   ";
		outputBobNames(stream, cluster.colorIndices, indexMap, site);
		stream << '\n';
	}
	stream << "Halfbrites\n";
	for (const Cluster &cluster : clusters) {
		if (!cluster.halfbriteIndices.empty()) {
			stream << counter++ << " - " << cluster.redMean/2 << ", " << cluster.greenMean/2 << ", " << cluster.blueMean/2 << "   -   ";
			outputBobNames(stream, cluster.halfbriteIndices, indexMap, site);
			stream << '\n';
		}
		else {
			counter++;
		}
	}
	stream << flush;
}

void extractColor(png_color &color, const Cluster &cluster, double multiplier=1.0) {
	color.red = min(255, static_cast<int>(cluster.redMean * multiplier));
	color.green = min(255, static_cast<int>(cluster.greenMean * multiplier));
	color.blue = min(255, static_cast<int>(cluster.blueMean * multiplier));
}

void extractPalette(png_colorp palette, const vector<Cluster> &clusters) {
	int counter = 0;
	for (const Cluster &cluster : clusters) {
		png_color &color = palette[counter++];
		extractColor(color, cluster);
	}
	for (const Cluster &cluster : clusters) {
		if (!cluster.halfbriteIndices.empty()) {
			png_color &color = palette[counter++];
			extractColor(color, cluster, 0.5);
		}
		else {
			counter++;
		}
	}
}

template<int cellW, int cellH, int nCellsW, int nCellsH>
struct DrawingContext {
	png_color rows[nCellsH][cellW*nCellsW] = { 0 };
	png_bytep pngRows[cellH*nCellsH];
	int cursX = 0, cursY = 0;
	int width = cellW*nCellsW;
	int height = cellH*nCellsH;
	DrawingContext() {
		for (size_t i=0; i < cellH * nCellsH; ++i) {
			pngRows[i] = (png_bytep)&(rows[i/cellH][0]);
		}
	}

	void tab() {
		cursX = nCellsW / 2;
	}

	void nextLine() {
		cursX = 0;
		cursY++;
	}

	void putColor(double red, double green, double blue) {
		if (cursX < 0 || nCellsW <= cursX) return;
		if (cursY < 0 || nCellsH <= cursY) return;
		png_color color;
		color.red = min(255, static_cast<int>(red));
		color.green = min(255, static_cast<int>(green));
		color.blue = min(255, static_cast<int>(blue));
		int columnStart = cursX*cellW;
		for (int columnOffset = 0; columnOffset < cellW; ++columnOffset) {
			rows[cursY][columnStart + columnOffset] = color;
		}
		cursX++;
	}

	void putColor(int color) {
		int red, green, blue;
		extractColorComponents(color, red, green, blue);
		putColor(red, green, blue);
	}
};

template <typename DrawingContext>
void drawClusterAverage(DrawingContext &canvas, const vector<Cluster> &clusters, const vector<int> colors) {
	for (const Cluster &cluster : clusters) {
		canvas.putColor(cluster.redMean, cluster.greenMean, cluster.blueMean);
		for (size_t index : cluster.colorIndices) {
			canvas.putColor(colors[index]);
		}
		canvas.tab();
		canvas.putColor(cluster.redMean / 2, cluster.greenMean / 2, cluster.blueMean / 2);
		for (size_t index : cluster.halfbriteIndices) {
			canvas.putColor(colors[index]);
		}
		canvas.nextLine();
	}
}

int pngClusterAverage(FILE *fp, const vector<Cluster> &clusters, const vector<int> colors) {
	png_structp writer = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!writer) {
		return errno;
	}
	png_infop info = png_create_info_struct(writer);
	if (!info) {
		png_destroy_write_struct(&writer, NULL);
		return errno;
	}

	if (setjmp(png_jmpbuf(writer))) {
		png_destroy_write_struct(&writer, &info);
		return errno;
	}

	png_init_io(writer, fp);

	DrawingContext<16, 16, 32, 32> canvas;
	drawClusterAverage(canvas, clusters, colors);

	png_set_IHDR(writer, info, canvas.width, canvas.height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_set_rows(writer, info, canvas.pngRows);

	png_write_png(writer, info, PNG_TRANSFORM_IDENTITY, NULL);

	return 0;
}

int pngClusters(FILE *fp, const vector<Cluster> &clusters) {
	png_structp writer = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!writer) {
		return errno;
	}
	png_infop info = png_create_info_struct(writer);
	if (!info) {
		png_destroy_write_struct(&writer, NULL);
		return errno;
	}

	if (setjmp(png_jmpbuf(writer))) {
		png_destroy_write_struct(&writer, &info);
		return errno;
	}

	png_init_io(writer, fp);

	png_byte row[16] = { 0 };
	png_bytep rows[16];
	for (png_bytep &initRow : rows) {
		initRow = row;
	}

	png_color palette[64] = { 0 };
	extractPalette(palette, clusters);

	png_set_IHDR(writer, info, 16, 16, 8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_set_PLTE(writer, info, palette, 64);
	png_set_rows(writer, info, rows);

	png_write_png(writer, info, PNG_TRANSFORM_IDENTITY, NULL);

	return 0;
}

void stiansShittyAlgorithm(vector<Cluster> &clusters, const vector<int> &allColors) {
	epsilonClusters(clusters, allColors, 1000.0);
}

void crudeTestAlgorithm(ostream &logfile, const Site &site) {
	const vector<int> &colors = site.colors;

	/*
	vector<int> allColors;
	for (int &color : colors) {
		allColors.push_back(color);
		allColors.push_back(halveColor(color));
	}
	vector<vector<int> > distances(allColors.size());
	vector<ColorDist> sortedDistances;
	for (size_t row = 0; row < allColors.size(); ++row) {
		for (size_t col = row + 1; col < allColors.size(); ++col) {
			int dist = squaredEuclideanDistance(allColors[row], allColors[col]);
			distances[row].push_back(dist);

			ColorDist colorDist;
			colorDist.distance = dist;
			colorDist.colorA = allColors[row];
			colorDist.colorB = allColors[col];
			colorDist.indexA = row;
			colorDist.indexB = col;
			sortedDistances.push_back(colorDist);
		}
	}

	sort(sortedDistances.begin(), sortedDistances.end());

	outputDistances(logfile, allColors, distances);
	logfile << '\n';
	outputLeastDistances(logfile, sortedDistances, 0, 10, true);
	logfile << '\n';
	outputLeastDistances(logfile, sortedDistances, 0, 10);
	logfile << '\n';
	outputLeastDistances(logfile, sortedDistances);

	*/

	vector<size_t> indexMap(colors.size());
	iota(indexMap.begin(), indexMap.end(), 0);

	sort(indexMap.rbegin(), indexMap.rend(), [&colors](size_t lhs, size_t rhs) {
		return colors[lhs] < colors[rhs];
	});

	vector<int> colorsSorted(colors.size());
	for (size_t i = 0; i < colorsSorted.size(); ++i) {
		colorsSorted[i] = colors[indexMap[i]];
	}

	vector<Cluster> clusters;

	stiansShittyAlgorithm(clusters, colorsSorted);
	outputClusters(logfile, clusters, indexMap, site);

	string pngfilename = site.name + ".png";
	FILE *fp;
	fopen_s(&fp, pngfilename.c_str(), "wb");
	errno = 0;
	int error = pngClusters(fp, clusters);
	fclose(fp);
	if (error) {
		char buf[1024];
		strerror_s(buf, error);
		perror(buf);
	}

	pngfilename = site.name + "Average.png";
	fopen_s(&fp, pngfilename.c_str(), "wb");
	errno = 0;
	error = pngClusterAverage(fp, clusters, colorsSorted);
	fclose(fp);
	if (error) {
		char buf[1024];
		strerror_s(buf, error);
		perror(buf);
	}
}

int main() {
	ifstream file("clarencecolors.csv");
	ofstream logfile("clarencelog.txt");
	vector<Site> sites;
	vector<Bob> bobs;
	bobs.reserve(200); // to avoid reallocation, since we will make pointers to elements in this vector.
	readfile(file, sites, bobs);

	crudeTestAlgorithm(logfile, sites[0]);
}