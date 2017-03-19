#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
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
			bob.name = cell + " - bg";
			bobs.push_back(bob);
			sites.back().bobs.push_back(&bobs.back());
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
			}
			// otherwise the cell contains the name of a bob
			else {

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
	for (int row = 0; row < distances.size(); ++row) {
		stream << setw(6) << hex << allColors[row] << " ";
		for (int col = 0; col <= row; ++col) {
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

struct Cluster {
	vector<size_t> colorIndeces;
	double redMean = 0.0;
	double greenMean = 0.0;
	double blueMean = 0.0;
	bool halfbrite = false;
	bool allHalfbrite = true;
	void insertColor(int color, size_t index) {
		size_t n = colorIndeces.size();
		int red = color >> 16;
		int green = color >> 8 & 0xFF;
		int blue = color & 0xFF;

		redMean = (redMean * n + red) / (n + 1);
		greenMean = (greenMean * n + green) / (n + 1);
		blueMean = (blueMean*n + blue) / (n + 1);
		colorIndeces.push_back(index);

		if (isHalfbrite(index)) {
			halfbrite = true;
		}
		else {
			allHalfbrite = false;
		}
	}

	double distance(int color) {
		int red = color >> 16;
		int green = color >> 8 & 0xFF;
		int blue = color & 0xFF;

		double dRed = red - redMean;
		double dGreen = green - greenMean;
		double dBlue = blue - blueMean;

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
		}
		if (!foundCluster) {
			Cluster cluster;
			cluster.insertColor(allColors[index], index);
			result.push_back(cluster);
		}
	}
}

void outputClusters(ostream &stream, const vector <Cluster> &clusters, const vector<int> &allColors) {
	int counter = 0;
	stream << "Colors\n";
	for (const Cluster &cluster : clusters) {
		if (!cluster.halfbrite) {
			stream << counter++ << " - " << cluster.redMean << ", " << cluster.greenMean << ", " << cluster.blueMean << '\n';
		}
	}
	stream << "Halfbrites\n";
	for (const Cluster &cluster : clusters) {
		if (cluster.halfbrite && !cluster.allHalfbrite) {
			stream << counter++ << " - " << cluster.redMean << ", " << cluster.greenMean << ", " << cluster.blueMean << '\n';
		}
	}
	stream << flush;
}

void stiansShittyAlgorithm(ostream &logfile, const vector<int> &allColors) {
	vector<Cluster> clusters;
	epsilonClusters(clusters, allColors, 900.0);
	outputClusters(logfile, clusters, allColors);
}

void crudeTestAlgorithm(ostream &logfile, vector<int> &colors) {
	vector<int> allColors;
	for (int &color : colors) {
		allColors.push_back(color);
		allColors.push_back(halveColor(color));
	}
	vector<vector<int> > distances(allColors.size());
	vector<ColorDist> sortedDistances;
	for (size_t row = 0; row < allColors.size(); ++row) {
		for (size_t col = row+1; col < allColors.size(); ++col) {
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


	stiansShittyAlgorithm(logfile, allColors);
}

int main() {
	ifstream file("clarencecolors.csv");
	ofstream logfile("clarencelog.txt");
	vector<Site> sites;
	vector<Bob> bobs;
	bobs.reserve(200); // to avoid reallocation, since we will make pointers to elements in this vector.
	readfile(file, sites, bobs);

	crudeTestAlgorithm(logfile, sites[0].colors);
}