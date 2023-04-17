//Author Raghav Gupta and Vikranth Gaddam
//CPP of $1 Recognizer
#include <wx/wx.h>
#include <wx/graphics.h>
#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <utility>
#include <stdlib.h>
#include <chrono>
#include <wx/math.h>
#include <iostream>
#include <tinyxml2.h>
#include <iostream>
#include <wx/log.h>
#include <vector>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <memory>
#include <unordered_map>
#include <wx/log.h>
#include <numeric>
#include <filesystem>
using namespace std;

class Point {//Defining point struct 
public:
	double x;
	double y;
	long long int t;
	Point(double x, double y) {
		this->x = x;
		this->y = y;
	}
	Point(double x, double y, int t) {
		this->x = x;
		this->y = y;
		this->t = t;
	}
	Point() {}
};


class splitDataStruct {
public:
	string name;
	vector<vector<Point>> gesturePoints; 
};

struct Rect { //Defining rectangle for bounding box
	double X;
	double Y;
	double Width;
	double Height;

	Rect(double x, double y, double width, double height)
		: X(x), Y(y), Width(width), Height(height) {}
};

const int NumMultistrokes = 16;
const int NumPoints = 96;
const double SquareSize = 250.0;
const double OneDThreshold = 0.25; // customize to desired gesture set (usually 0.20 - 0.35)
const Point Origin(0, 0);
const double Diagonal = sqrt(SquareSize * SquareSize + SquareSize * SquareSize);
const double HalfDiagonal = 0.5 * Diagonal;
const double AngleRange = 45.0 * M_PI / 180.0;
const double AnglePrecision = 2.0 * M_PI / 180.0;
const double Phi = 0.5 * (-1.0 + sqrt(5.0)); // Golden Ratio
const int StartAngleIndex = (NumPoints / 8); // eighth of gesture length
const double AngleSimilarityThreshold = 30.0 * M_PI / 180.0;

//double Deg2Rad(double d) {
//	return (d * M_PI / 180.0);
//}



struct Unistroke {//This is where preprocessing happens for points
public:
	std::string name;
	std::vector<Point> points;
	std::vector<double> vectorizedPoints;
	Point StartUnitVector;
	Unistroke() {}
	Unistroke(string& name, bool useBoundedRotationInvariance, vector<Point>& points) //Constructor where the point are preprocessed and stored
	{
		//Note 1.Resample 2.Rotate By needs indicative angle 3.Scale 4.Transalate 5.For Protractor -- Vectorize
		this->name = name;
		this->points = Resample(points, NumPoints);//1
		double radians = IndicativeAngle(this->points);//2
		this->points = RotateBy(this->points, -radians);//3
		this->points = ScaleDimTo(this->points, SquareSize, OneDThreshold);//4
		if (useBoundedRotationInvariance) {
			points = RotateBy(this->points, radians);
		}
		this->points = TranslateTo(this->points, Origin);//5
		this->StartUnitVector = CalcStartUnitVector(this->points, StartAngleIndex); //should check it again
		//wxLogMessage("%f ,%f",StartUnitVector.x,StartUnitVector.y);
		this->vectorizedPoints = Vectorize(this->points); // for Protractor	
	}
	std::vector<Point>Resample(std::vector<Point>& points, int n)
	{

		double I = PathLength(points) / (n - 1); // interval length
		double D = 0.0;
		vector<Point> newpoints = { points[0] };
		for (int i = 1; i < points.size(); i++)
		{
			double d = Distance(points[i - 1], points[i]);
			if ((D + d) >= I)
			{
				double qx = points[i - 1].x + ((I - D) / d) * (points[i].x - points[i - 1].x);
				double qy = points[i - 1].y + ((I - D) / d) * (points[i].y - points[i - 1].y);
				Point q = Point(qx, qy);
				newpoints.push_back(q); // append new point 'q'
				points.insert(points.begin() + i, q); // insert 'q' at position i in points s.t. 'q' will be the next i
				D = 0.0;
			}
			else D += d;
		}
		if (newpoints.size() == n - 1) // somtimes we fall a rounding-error short of adding the last point, so add it if so
			newpoints.push_back(Point(points[points.size() - 1].x, points[points.size() - 1].y));
		return newpoints;
		//TODO Vikranth
	}
	double IndicativeAngle(std::vector<Point> points)
	{
		Point c = Centroid(points);
		return atan2(c.y - points[0].y, c.x - points[0].x);
	}

	vector<Point> ScaleDimTo(vector<Point> points, double size, double ratio1D)
	{
		Rect B = BoundingBox(points);
		bool uniformly = min(B.Width / B.Height, B.Height / B.Width) <= ratio1D;
		vector<Point> newpoints;
		for (int i = 0; i < points.size(); i++) {
			double qx = uniformly ? points[i].x * (size / max(B.Width, B.Height)) : points[i].x * (size / B.Width);
			double qy = uniformly ? points[i].y * (size / max(B.Width, B.Height)) : points[i].y * (size / B.Height);
			newpoints.push_back(Point(qx, qy));
		}
		return newpoints;
	}

	std::vector<Point> TranslateTo(const std::vector<Point>& points, Point pt)
	{
		Point c = Centroid(points);
		std::vector<Point> newPoints;
		for (const auto& point : points)
		{
			Point q = { point.x + pt.x - c.x, point.y + pt.y - c.y };
			newPoints.push_back(q);
		}//TODO Vikranth
		return newPoints;

	}
	std::vector<double> Vectorize(const std::vector<Point>& points) {
		double sum = 0.0;
		std::vector<double> vector;
		for (const auto& point : points) {
			vector.push_back(point.x);
			vector.push_back(point.y);
			sum += point.x * point.x + point.y * point.y;
		}//TODO Raghav
		double magnitude = sqrt(sum);
		for (auto& val : vector)
			val /= magnitude;
		return vector;
	}

	std::vector<Point> RotateBy(std::vector<Point>& points, double radians) // rotates points around centroid
	{//TODO Vikranth
		Point c = Centroid(points);
		double cos = std::cos(radians);
		double sin = std::sin(radians);
		std::vector<Point> newpoints;
		for (auto& point : points) {
			double qx = (point.x - c.x) * cos - (point.y - c.y) * sin + c.x;
			double qy = (point.x - c.x) * sin + (point.y - c.y) * cos + c.y;
			newpoints.emplace_back(qx, qy);
		}
		return newpoints;
	}


	Rect BoundingBox(std::vector<Point> points)
	{//TODO Vikranth
		double minX = INFINITY, maxX = -INFINITY, minY = INFINITY, maxY = -INFINITY;
		for (int i = 0; i < points.size(); i++) {
			minX = (double)std::min(minX, points[i].x);
			minY = (double)std::min(minY, points[i].y);
			maxX = (double)std::max(maxX, points[i].x);
			maxY = (double)std::max(maxY, points[i].y);
		}
		return Rect(minX, minY, maxX - minX, maxY - minY);
	}
	double PathLength(std::vector<Point> points)
	{
		double d = 0.0;
		for (int i = 1; i < points.size(); i++)
			d += Distance(points[i - 1], points[i]);
		return d;
	}
	double Distance(Point p1, Point p2)
	{
		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;
		return sqrt(dx * dx + dy * dy);
	}
	Point Centroid(const std::vector<Point>& points)
	{
		double x = 0.0, y = 0.0;
		for (const auto& point : points) {
			x += point.x;
			y += point.y;
		}
		x /= points.size();
		y /= points.size();
		return { x, y };
	}

	Point CalcStartUnitVector(vector<Point> points, int index) {
		Point v(points[index].x - points[0].x, points[index].y - points[0].y);
		double len = sqrt(v.x * v.x + v.y * v.y);
		return Point(v.x / len, v.y / len);
	}



};


struct Multistroke {//This is where preprocessing happens for points
public:
	std::string name;
	std::vector<Point> points;
	std::vector<double> vectorizedPoints;
	int NumStrokes;
	std::vector<Unistroke> Unistrokes;
	Multistroke() {}

	Multistroke(std::string name, bool useBoundedRotationInvariance, std::vector<std::vector<Point>> strokes) {
		this->name = name;
		this->NumStrokes = strokes.size();

		vector<int> order(strokes.size());
		for (int i = 0; i < strokes.size(); i++) {
			order[i] = i;
		}

		std::vector<std::vector<int>> orders;
		HeapPermute(strokes.size(), order, orders);

		std::vector<std::vector<Point>> unistrokes = MakeUnistrokes(strokes, orders);
		this->Unistrokes.resize(unistrokes.size());

		for (int j = 0; j < unistrokes.size(); j++) {
			this->Unistrokes[j] = Unistroke(name, useBoundedRotationInvariance, unistrokes[j]);
		}
	}

private:
	void HeapPermute(int n, vector<int>& order, vector<vector<int>>& orders) {
		if (n == 1) {
			orders.push_back(order);
		}
		else {
			for (int i = 0; i < n; i++) {
				HeapPermute(n - 1, order, orders);
				if (n % 2 == 1) { // swap 0, n-1
					int tmp = order[0];
					order[0] = order[n - 1];
					order[n - 1] = tmp;
				}
				else { // swap i, n-1
					int tmp = order[i];
					order[i] = order[n - 1];
					order[n - 1] = tmp;
				}
			}
		}
	}
	vector<vector<Point>> MakeUnistrokes(vector<vector<Point>> strokes, vector<vector<int>> orders)
	{
		vector<vector<Point>> unistrokes; // array of point arrays
		for (int r = 0; r < orders.size(); r++)
		{
			for (int b = 0; b < pow(2, orders[r].size()); b++) // use b's bits for directions
			{
				vector<Point> unistroke; // array of points
				for (int i = 0; i < orders[r].size(); i++)
				{
					vector<Point> pts;
					if (((b >> i) & 1) == 1) // is b's bit at index i on?
						pts = vector<Point>(strokes[orders[r][i]].rbegin(), strokes[orders[r][i]].rend()); // copy and reverse
					else
						pts = strokes[orders[r][i]]; // copy
					for (int p = 0; p < pts.size(); p++)
						unistroke.push_back(pts[p]); // append points
				}
				unistrokes.push_back(unistroke); // add one unistroke to set
			}
		}
		return unistrokes;
	}
};


struct Result {//Result struct for displaying in canvas
	std::string Name;
	double Score;
	int Time;

	Result(const std::string& name, double score, int time) : Name(name), Score(score), Time(time) {}
};

struct OfflineResult {
	string gestureName;
	double score = 0.0;
	vector<pair<string, double>> nbest;
};

struct LogData {
	string User;
	string GestureType;
	string RandomIteration;
	int NoTrainingExample;
	int TotalSizeOfTrainingSet;
	vector<string>TrainingSetContents;
	string candidateSpecificInstance;
	string RecoResultGestureType;//what was recognized
	bool correct;
	string RecoResultScore;
	string RecoResultBestMatch;
	vector<pair<string, double>> RecoResultNBest;
};


class MGestureRecognizer {

	
public:
	MGestureRecognizer() {}
	MGestureRecognizer(bool useBoundedRotationInvariance)
	{
		// one predefined multistroke for each multistroke type
		vector<vector<Point>> Tpoints = {
			{Point(30, 7), Point(103, 7)},
		{Point(66, 7), Point(66, 87)}
		};
		vector<vector<Point>> Npoints = {
			{Point(177,92),Point(177,2)},
		{Point(182,1), Point(246,95)},
		{Point(247,87), Point(247,1)}
		};
		vector<vector<Point>> Dpoints = {
		{Point(345,9),Point(345,87)},
		{Point(351,8),Point(363,8),Point(372,9),Point(380,11),Point(386,14),Point(391,17),Point(394,22),Point(397,28),Point(399,34),Point(400,42),Point(400,50),Point(400,56),Point(399,61),Point(397,66),Point(394,70),Point(391,74),Point(386,78),Point(382,81),Point(377,83),Point(372,85),Point(367,87),Point(360,87),Point(355,88),Point(349,87)}
		};
		vector<vector<Point>> Ppoints = {
		{Point(507,8),Point(507,87)},
		{Point(513,7),Point(528,7),Point(537,8),Point(544,10),Point(550,12),Point(555,15),Point(558,18),Point(560,22),Point(561,27),Point(562,33),Point(561,37),Point(559,42),Point(556,45),Point(550,48),Point(544,51),Point(538,53),Point(532,54),Point(525,55),Point(519,55),Point(513,55),Point(510,55)}
		};

		vector<vector<Point>> Xpoints = {
			{Point(30,146),Point(106,222)},
			{Point(30,225),Point(106,146)}
		};
		vector<vector<Point>> Hpoints = {
			{Point(188,137),Point(188,225)},
			{Point(188,180),Point(241,180)},
			{Point(241,137),Point(241,225)}
		};

		vector<vector<Point>> Ipoints = {
			{Point(371,149),Point(371,221)},
			{Point(341,149),Point(401,149)},
			{Point(341,221),Point(401,221)}
		};

		vector<vector<Point>> exclamationspoints = {
			{Point(526,142),Point(526,204)},
			{Point(526,221)}
		};
		vector<vector<Point>> linePoints = {
		{Point(12,347), Point(119,347)}
		};

		vector<vector<Point>> fivePointStarPoints = {
			{Point(177,396), Point(223,299), Point(262,396), Point(168,332), Point(278,332), Point(184,397)}
		};

		vector<vector<Point>> nullPoints = {
			{Point(382,310),Point(377,308),Point(373,307),Point(366,307),Point(360,310),Point(356,313),Point(353,316),Point(349,321),Point(347,326),Point(344,331),Point(342,337),Point(341,343),Point(341,350),Point(341,358),Point(342,362),Point(344,366),Point(347,370),Point(351,374),Point(356,379),Point(361,382),Point(368,385),Point(374,387),Point(381,387),Point(390,387),Point(397,385),Point(404,382),Point(408,378),Point(412,373),Point(416,367),Point(418,361),Point(419,353),Point(418,346),Point(417,341),Point(416,336),Point(413,331),Point(410,326),Point(404,320),Point(400,317),Point(393,313),Point(392,312)},
			{Point(418,309), Point(337,390)}
		};

		vector<vector<Point>> arrowheadPoints = {
			{Point(506,349), Point(574,349)},
			{Point(525,306), Point(584,349), Point(525,388)}
		};

		vector<vector<Point>> pitchforkPoints = {
			{Point(38,470), Point(36,476), Point(36,482), Point(37,489), Point(39,496), Point(42,500), Point(46,503), Point(50,507), Point(56,509), Point(63,509), Point(70,508), Point(75,506), Point(79,503), Point(82,499), Point(85,493), Point(87,487), Point(88,480), Point(88,474), Point(87,468)},
			{Point(62,464), Point(62,571)}
		};

		vector<vector<Point>> sixPointStarPoints = {
			{Point(177,554), Point(223,476), Point(268,554), Point(183,554)},
			{Point(177,490), Point(223,568), Point(268,490), Point(183,490)}
		};

		vector<vector<Point>> asteriskPoints = {
			{Point(325,499), Point(417,557)},
			{Point(417,499), Point(325,557)},
			{Point(371,486), Point(371,571)}
		};
		vector<vector<Point>> half_notepoints = {
	{Point(546, 465), Point(546, 531)},
	{Point(540, 530), Point(536, 529), Point(533, 528), Point(529, 529), Point(524, 530), Point(520, 532), Point(515, 535), Point(511, 539), Point(508, 545), Point(506, 548), Point(506, 554), Point(509, 558), Point(512, 561), Point(517, 564), Point(521, 564), Point(527, 563), Point(531, 560), Point(535, 557), Point(538, 553), Point(542, 548), Point(544, 544), Point(546, 540), Point(546, 536)}
		};

		vector<string> gestures = { "T", "N", "D", "P", "X", "H", "I", "exclamation"
							  "line", "five-point star", "null", "arrowhead", "pitchfork","six-point star","asterisk", "half-note" };

		Multistrokes.push_back(Multistroke("T", useBoundedRotationInvariance, { Tpoints }));
		 Multistrokes.push_back(Multistroke("N", useBoundedRotationInvariance, { Npoints }));
		 Multistrokes.push_back(Multistroke("D", useBoundedRotationInvariance, { Dpoints }));
		 Multistrokes.push_back(Multistroke("P", useBoundedRotationInvariance, { Ppoints }));
		 Multistrokes.push_back(Multistroke("X", useBoundedRotationInvariance, { Xpoints }));
		 Multistrokes.push_back(Multistroke("H", useBoundedRotationInvariance, { Hpoints }));
		 Multistrokes.push_back(Multistroke("I", useBoundedRotationInvariance, { Ipoints }));
		Multistrokes.push_back(Multistroke("exclamation", useBoundedRotationInvariance, { exclamationspoints }));
		 Multistrokes.push_back(Multistroke("line", useBoundedRotationInvariance, { linePoints }));
		 Multistrokes.push_back(Multistroke("five-point star", useBoundedRotationInvariance, { fivePointStarPoints }));
		Multistrokes.push_back(Multistroke("null", useBoundedRotationInvariance, { nullPoints }));
		Multistrokes.push_back(Multistroke("arrowhead", useBoundedRotationInvariance, { arrowheadPoints }));
		 Multistrokes.push_back(Multistroke("pitchfork", useBoundedRotationInvariance, { pitchforkPoints }));
		 Multistrokes.push_back(Multistroke("six-point star", useBoundedRotationInvariance, { sixPointStarPoints }));
		 Multistrokes.push_back(Multistroke("asterisk", useBoundedRotationInvariance, { asteriskPoints }));
		 Multistrokes.push_back(Multistroke("half-note", useBoundedRotationInvariance, { half_notepoints }));

	}
	vector<Multistroke> OfflineStrokes;
	MGestureRecognizer(vector<splitDataStruct>& inputTempatePoints) {
		
		for (auto& it : inputTempatePoints) {
			//wxLogMessage("Preprcoess started for %s", it.name);
			this->OfflineStrokes.push_back(Multistroke(it.name,false,it.gesturePoints));
			//wxLogMessage("Preprcoess Done for %s",it.name);
		}
		//wxLogMessage("Preprcoess DOne");
	}
	bool comparePair(const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
		return a.second > b.second;
	}

	vector<Multistroke> Multistrokes;
	Result Recognize(std::vector<std::vector<Point>> strokes, bool useBoundedRotationInvariance, bool requireSameNoOfStrokes, bool useProtractor) {
		auto t0 = std::chrono::high_resolution_clock::now();
		auto points = CombineStrokes(strokes); // make one connected unistroke from the given strokes
		string temp = "";
		Unistroke candidate(temp, useBoundedRotationInvariance, points);

		int u = -1;
		int i = 0;
		int j = 0;
		double b = INFINITY;
		for (i = 0; i < Multistrokes.size(); i++) // for each multistroke template
		{
			if (!requireSameNoOfStrokes || strokes.size() == Multistrokes[i].NumStrokes) // optional -- only attempt match when same # of component strokes
			{
				for (j = 0; j < Multistrokes[i].Unistrokes.size(); j++) // for each unistroke within this multistroke
				{
					if (AngleBetweenUnitVectors(candidate.StartUnitVector, Multistrokes[i].Unistrokes[j].StartUnitVector) <= AngleSimilarityThreshold) // strokes start in the same direction
					{
						double d;
						if (useProtractor)
							d = OptimalCosineDistance(Multistrokes[i].Unistrokes[j].vectorizedPoints, candidate.vectorizedPoints);
						else
							d = DistanceAtBestAngle(candidate.points, Multistrokes[i].Unistrokes[j], -AngleRange, AngleRange, AnglePrecision);
						if (d < b) {
							b = d;
							u = i;
						}
					}
				}
			}
		}
		auto t1 = std::chrono::high_resolution_clock::now();
		return (u == -1) ? Result("No match.", 0.0, std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) : Result(Multistrokes[u].name, useProtractor ? (1.0 - b) : (1.0 - b / HalfDiagonal), std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
	};
	OfflineResult OfflineRecognize(vector<vector<Point>> strokes, bool useBoundedRotationInvariance, bool requireSameNoOfStrokes, bool useProtractor) {
		//wxLogMessage("Here\n");
		auto t0 = std::chrono::high_resolution_clock::now();
		/*wxLogMessage("Here after\n");*/
		vector<OfflineResult> final;
		auto points = CombineStrokes(strokes); // make one connected unistroke from the given strokes
		string t = "";
		//wxLogMessage("Reading candidatae");
		Unistroke candidate(t, useBoundedRotationInvariance, points);

		int u = -1;
		int i = 0;
		int j = 0;
		double b = INFINITY;
		OfflineResult temp;
		vector<pair<string, double>> t1;
		for (i = 0; i < OfflineStrokes.size(); i++) // for each multistroke template
		{
			if (!requireSameNoOfStrokes || strokes.size() == OfflineStrokes[i].NumStrokes) // optional -- only attempt match when same # of component strokes
			{
				for (j = 0; j < OfflineStrokes[i].Unistrokes.size(); j++) // for each unistroke within this multistroke
				{
					if (AngleBetweenUnitVectors(candidate.StartUnitVector, OfflineStrokes[i].Unistrokes[j].StartUnitVector) <= AngleSimilarityThreshold) // strokes start in the same direction
					{
						
						double d=0.0;
						if (useProtractor)
							d = OptimalCosineDistance(OfflineStrokes[i].Unistrokes[j].vectorizedPoints, candidate.vectorizedPoints);
						else
						{
							d = DistanceAtBestAngle(candidate.points, OfflineStrokes[i].Unistrokes[j], -AngleRange, AngleRange, AnglePrecision);
						}
							if (d < b) {
								b = d;
								u = i;
							}
							temp.gestureName = OfflineStrokes[i].name;
							temp.score = 1.0 - d / HalfDiagonal;
							string name = temp.gestureName;
							t1.push_back(make_pair(name, temp.score));
					}
				}
			}
		}

		if (u != -1) {
			temp.gestureName = OfflineStrokes[u].name;
		}
		else {
			temp.gestureName = "No match.";
		}
		temp.score = (1.0 - b / HalfDiagonal);
		std::sort(t1.begin(), t1.end(), std::bind(&MGestureRecognizer::comparePair, this, std::placeholders::_1, std::placeholders::_2));
		const int max_size = 50;
		std::vector<std::pair<std::string, double>> t2(max_size);

		// copy the first 50 elements from t1 to t2
		const int size = std::min(static_cast<int>(t1.size()), max_size);
		std::copy_n(t1.begin(), size, t2.begin());
		temp.nbest = t2;
		t2.clear();
		t1.clear();
		//auto t1 = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();//Time elapsed
		return temp;


		/*auto t1 = std::chrono::high_resolution_clock::now();
		return (u == -1) ? Result("No match.", 0.0, std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) : Result(OfflineStrokes[u].name, useProtractor ? (1.0 - b) : (1.0 - b / HalfDiagonal), std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());*/

	}


	double OptimalCosineDistance(const std::vector<double>& v1, const std::vector<double>& v2)
	{
		double a = 0.0;
		double b = 0.0;

		for (int i = 0; i < v1.size(); i += 2)
		{
			a += v1[i] * v2[i] + v1[i + 1] * v2[i + 1];
			b += v1[i] * v2[i + 1] - v1[i + 1] * v2[i];
		}

		double angle = atan(b / a);
		return acos(a * cos(angle) + b * sin(angle));
	}
	double DistanceAtBestAngle(std::vector<Point>& points, Unistroke& T, double a, double b, double threshold) {
		//TODO Vikranth
		double x1 = Phi * a + (1.0 - Phi) * b;
		double f1 = DistanceAtAngle(points, T, x1);
		double x2 = (1.0 - Phi) * a + Phi * b;
		double f2 = DistanceAtAngle(points, T, x2);
		while (std::abs(b - a) > threshold) {
			if (f1 < f2) {
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = Phi * a + (1.0 - Phi) * b;
				f1 = DistanceAtAngle(points, T, x1);
			}
			else {
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = (1.0 - Phi) * a + Phi * b;
				f2 = DistanceAtAngle(points, T, x2);
			}
		}
		return std::min(f1, f2);
	}
	double DistanceAtAngle(std::vector<Point>& points, Unistroke& T, double radians)
	{
		//TODO Raghav
		std::vector<Point> newpoints = RotateBy(points, radians);

		return PathDistance(newpoints, T);
	}
	std::vector<Point> RotateBy(std::vector<Point>& points, double radians) // rotates points around centroid
	{
		//TODO Vikranth
		Point c = Centroid(points);
		double cos = std::cos(radians);
		double sin = std::sin(radians);
		std::vector<Point> newpoints;
		for (auto& point : points) {
			double qx = (point.x - c.x) * cos - (point.y - c.y) * sin + c.x;
			double qy = (point.x - c.x) * sin + (point.y - c.y) * cos + c.y;
			newpoints.emplace_back(qx, qy);
		}
		return newpoints;
	}

	double PathDistance(vector<Point> pts1, Unistroke inputStrokes)

	{//TODO Vikranth
		std::vector<Point> pts2 = inputStrokes.points;//TODO 
		double d = 0.0;
		for (int i = 0; i < pts1.size(); i++) // assumes pts1.size() == pts2.size()
			d += Distance(pts1[i], pts2[i]);
		return d / pts1.size();
	}

	Point Centroid(const std::vector<Point>& points)
	{
		//TODO Vikranth
		double x = 0.0, y = 0.0;
		for (const auto& point : points) {
			x += point.x;
			y += point.y;
		}
		x /= points.size();
		y /= points.size();
		return { x, y };
	}
	double Distance(Point p1, Point p2)
	{
		//TODO Vikranth
		double dx = p2.x - p1.x;
		double dy = p2.y - p1.y;
		return std::sqrt(dx * dx + dy * dy);
	}
	double AngleBetweenUnitVectors(Point v1, Point v2) {
		double n = (v1.x * v2.x + v1.y * v2.y);
		double c = fmax(-1.0, fmin(1.0, n)); // ensure [-1,+1]
		return acos(c); // arc cosine of the vector dot product
	}
	int AddGesture(std::string name, bool useBoundedRotationInvariance, std::vector<std::vector<Point>> strokes) {
		Multistrokes.emplace_back(name, useBoundedRotationInvariance, strokes);
		int num = 0;
		for (int i = 0; i < Multistrokes.size(); i++) {
			if (Multistrokes[i].name == name)
				num++;
		}
		return num;
	}
	vector<Point> CombineStrokes(vector<vector<Point>> strokes) {
		vector<Point> points;
		for (int s = 0; s < strokes.size(); s++) {
			for (int p = 0; p < strokes[s].size(); p++) {
				points.push_back(Point(strokes[s][p].x, strokes[s][p].y));
			}
		}
		return points;
	}


	int DeleteUserGestures() {
		Multistrokes.resize(NumMultistrokes); // clear any beyond the original set
		return NumMultistrokes;
	}
};



class OfflineRecognizer {
public:
	//GestureRecognizer recognizer;
	map< string, map< string, map< string, vector< vector<Point>>>>> offlineData;
	//[s02][""medium]["arrow01][]
	map< string, map< string, map< string, vector< vector<Point>>>>> preProcessedData;
	//Gesture name  list of points map("arrow01",arrow01,points[]) user01 speed gestureName points

	//$1n algo
	//GestureRecognizer recognizer;
	

	//Input points from XML Files 
	OfflineRecognizer() {
		tinyxml2::XMLDocument doc;
		wxLog::SetActiveTarget(new wxLogStderr);
		//XML File 
		vector<string> labelList = { "arrowhead", "asterisk", "D", "exclamation_point", "five_point_star", "H", "half_note", "I"
										  ,"line", "N", "null", "P", "pitchfork","six_point_star","T", "X" };
		vector<string> userName = { "10" ,"11","12","22","28","41","58","61","66","68","71","73","75","77","85","88","94","95","98","99"};
		string fileName = "";
		string totalFileName = "";
		string part = "";
		string storeName = "";
		part = ""; //C:/Users/vikra/Desktop/HCIRA/xml/xml_logs/s0
		//C:\Users\vikra\Desktop\0703\ndollar\10-stylus-MEDIUM
		for (const std::string& name : userName) {
			for (int i = 0; i < 16; i++) {
				string temp = "";
				if (name=="94"||name == "66" || name == "58" || name == "61" || name == "73" || name == "75" || name == "77" || name == "99" || name == "41" || name == "85") {
					temp = "-finger-MEDIUM";
				}
				else {
					temp = "-stylus-MEDIUM";
				}
				part = "C:/Users/91977/Downloads/mmg(1)/" + name + temp + "/" + name + temp+ "-" + labelList[i] + "-";
				//wxLogMessage("Part :%s", part);
				for (int j = 1; j <= 10; j++) {
					vector<vector<Point>> strokes;
					if (j < 10) {
						fileName = part + "0" + to_string(j) + ".xml";
						totalFileName = labelList[i] + "0" + to_string(j);
					}
					else {
						fileName = part + to_string(j) + ".xml";
						totalFileName = labelList[i] + to_string(j);
					}

					doc.LoadFile(fileName.c_str());
					//wxLogMessage("File name %s, %s", fileName,totalFileName);

					tinyxml2::XMLElement* gesture = doc.FirstChildElement("Gesture");

					tinyxml2::XMLElement* stroke = gesture->FirstChildElement("Stroke");
					while (stroke != nullptr)
					{
						vector<Point> temp;
						tinyxml2::XMLElement* point = stroke->FirstChildElement("Point");
						while (point != nullptr)
						{
							double x = point->IntAttribute("X");
							double y = point->IntAttribute("Y");
							Point tempPoint(x, y);
							temp.push_back(tempPoint);
							point = point->NextSiblingElement("Point");
							
						}
						strokes.push_back(temp);
						stroke = stroke->NextSiblingElement("Stroke");
					}
					offlineData[name]["medium"][totalFileName] = strokes;
				}
				
			}
		}

	}
	//username finger medium arrow 01 vector<Vector<Point> for  single fiile 

	//void preProcessOfflineData() {
	//	preProcessedData = offlineData;

	//	for (auto& user : offlineData) {
	//		for (auto& speed : user.second) {
	//			if (speed.first == "medium") {
	//				for (auto& gesture : speed.second) {
	//					//wxLogMessage("");
	//					preProcessedData[user.first][speed.first][gesture.first] = vector< vector<Point>>();
	//					for (auto& temp : gesture.second) {
	//						//gesture name , points
	//						string s = gesture.first;
	//						Unistroke unistroke(s, false, temp);
	//						preProcessedData[user.first][speed.first][gesture.first].push_back({ unistroke.points });
	//					}
	//				}
	//			}
	//		}
	//	}
	//	for (auto& user : preProcessedData) {
	//		wxLogMessage("User: %s", user.first);
	//		for (auto& speed : user.second) {
	//			wxLogMessage("  Speed: %s", speed.first);
	//			for (auto& gesture : speed.second) {
	//				wxLogMessage("    Gesture: %s", gesture.first);
	//				//for (auto& points : gesture.second) {
	//				//	//wxLogMessage("      Points: ");
	//				//	for (auto& point : points) {
	//				//		//wxLogMessage("        x: %f, y: %f", point.x, point.y);
	//				//	}
	//				//}
	//			}
	//		}
	//	}

	//}

	void recognizeOfflineData() {

		/*for (auto& it : offlineData) {
			wxLogMessage("Preprcoess started for %s", it.first);
			for (auto& it2 : it.second) {
				wxLogMessage("speed %s", it2.first);
				for (auto& it3 : it2.second) {
					wxLogMessage("GEsture Name started for %s", it3.first);
					for (auto& it4 : it3.second) {
					
						for (auto& i : it4) {
							wxLogMessage("%f,%f", i.x, i.y);
						}
						}
				}
			}
			this->OfflineStrokes.push_back(Multistroke(it.name,false,it.gesturePoints));
			wxLogMessage("Preprcoess Done for %s", it.first);
		}*/

		//n best list sort of 
		/*wxLogMessage("Inside the recognizeOffline Data\n");
		wxLogMessage("Size map: %d", offlineData.size());*/
		map< string, map<int, map< string, double>>> score;//loading 
		double total = 0;
		double correct = 0;

		std::ofstream outfile("output.csv");
		LogData log;

		// Check if the file was opened successfully
		if (!outfile.is_open()) {
			std::cerr << "Failed to open file for writing." << std::endl;
		}

		//User[all-users]	GestureType[all-gestures-types]	RandomIteration[1to100]	#ofTrainingExamples[E]	TotalSizeOfTrainingSet[count]	TrainingSetContents[specific-gesture-instances]	Candidate[specific-instance]	RecoResultGestureType[what-was-recognized]	CorrectIncorrect[1or0]	RecoResultScore	RecoResultBestMatch[specific-instance]	RecoResultNBestSorted[instance-and-score]

		outfile << "User[all-users],GestureType[all-gesture-types],RandomIteration[1to10],#ofTrainingExample[E],TotalSizeOfTrainingSet[count],TrainingSetContents[specific-gesture-instances],Candidate[specific-instance],RecoResultGestureType[what-was-recognized],CorrectIncorrect[1or0],RecoResultScore,RecoResultBestMatch[specific-instance],RecoResultNBestSorted[instance-and-score]" << std::endl;
		for (auto& user : offlineData) {
			//wxLogMessage("Users data%s", user.first);
			//if (user.first == "s02")
			//wxLogMessage("users\n");
			{
				log.User = user.first;
				//score[user.first] = map<int,map< string, double>>();
				for (int gesture = 1; gesture <= 9; gesture++) { //this is the value of the 
					//wxLogMessage("gesture\n");
					//score[user.first][gesture] = map<string, double>();
					//score[user.first][gesture]["accuracy"] = 0.0;
					for (int i = 1; i <= 9; i++) { //we can decrease the value to 10
						//wxLogMessage("for loop\n");
						log.RandomIteration = to_string(i);
						log.NoTrainingExample = gesture;
						

						
						pair<splitDataStruct,vector<splitDataStruct> > split_data = getsplitdata(offlineData[user.first]["medium"], gesture);
						 vector<splitDataStruct> train_set = split_data.second;//TEsting conatines everything from offline 160-training
						splitDataStruct test_set = split_data.first;//a specific set of gestures with points training
					/*	wxLogMessage("Candidate name: %s",test_set.name	);
						wxLogMessage("training Size: %d", train_set.size());*/
						log.TotalSizeOfTrainingSet = train_set.size();
						log.TrainingSetContents.clear();
						for (auto& elem : train_set) {
							//user name // arrow 04 //iteration number
							string temp = user.first + "-" + elem.name + "-" + to_string(i);
							log.TrainingSetContents.push_back(temp);
						}

						
							//log.candidateSpecificInstance = user.first + "-" + test_set.name + "-" + to_string(i);
						//wxLogMessage("train set size: %d", train_set.size());
					

						MGestureRecognizer recognizer(train_set);//loading templates
						std::string r1 = ""; //gesturType

						for (char c : test_set.name) {
							if (!isdigit(c)) {
								r1 += c;
							}
						}

							
					OfflineResult res = recognizer.OfflineRecognize(test_set.gesturePoints, true, false, false);

							
							//wxLogMessage("Result: %s -  score: %f", res.Name, res.Score);
						string gestureRecognized = res.gestureName;
						vector<pair<string, double>> Nbest = res.nbest;

						std::string r2 = ""; //Recoresult gesture type

						for (char c : gestureRecognized) {
							if (!isdigit(c)) {
								r2 += c;
							}
						}
						log.GestureType = r1;
						log.RecoResultGestureType = r2;
						if (log.GestureType == log.RecoResultGestureType) {
							log.correct = true;
						}
						else {
							log.correct = false;
						}
						//log.correct = gestureRecognized.substr(0, elem.first.length() - 1) == elem.first.substr(0, elem.first.length() - 1);
						log.RecoResultScore = to_string(res.score);
						log.RecoResultBestMatch = user.first + "-" + res.gestureName;
						log.RecoResultNBest = res.nbest;
						log.candidateSpecificInstance = user.first + "-" + test_set.name;

						if (log.GestureType == log.RecoResultGestureType) {

							score[user.first][gesture][test_set.name] += 1;

							//log.GestureType = elem.first;
							correct += 1.0;

						}

						total += 1.0;


						std::string contentsString = vector1ToString(log.TrainingSetContents);

						std::string nbestStr = vectorToString(log.RecoResultNBest);
						outfile << log.User << "," << log.GestureType << "," << log.RandomIteration << "," << log.NoTrainingExample << "," << log.TotalSizeOfTrainingSet << "," << contentsString << "," << log.candidateSpecificInstance << "," << log.RecoResultGestureType << "," << to_string(log.correct) << "," << log.RecoResultScore << "," << log.RecoResultBestMatch << "," << nbestStr << endl;

						

					}

					//score[user.first][gesture]["accuracy"] =(correct/total) * 100.0;
				}
			}

		}
	/*	wxLogMessage("total %f", total);
		wxLogMessage("correct %f", correct);*/

		outfile << "The Total Accuracy" << endl;
		outfile << to_string((correct / total) * 100.0) << endl;
		/*for (const auto& user : score) {
			wxLogMessage("User: %s", user.first);
			for (const auto& gesture : user.second) {
				wxLogMessage("  Gesture: %d", gesture.first);
				for (const auto& data : gesture.second) {
					wxLogMessage("    Data: %s %f", data.first, data.second);
				}
			}
		}*/
		wxLogMessage("closing");
		outfile.close();
		//bool result=writeToFile(score);
		//wxLogMessage("Print successfully done %s",to_string(result));
		// calculate and output average accuracy

	}
	std::string vector1ToString(const std::vector<std::string>& v) {
		std::stringstream ss;
		for (const auto& s : v) {
			ss << s << ";";
		}
		return ss.str();
	}

	std::string vectorToString(const std::vector<std::pair<std::string, double>>& v) {
		std::stringstream ss;
		for (const auto& p : v) {
			ss << p.first << ";" << p.second << ";";
		}
		return ss.str();
	}

	pair<splitDataStruct,vector<splitDataStruct>> getsplitdata(map<string, vector<vector<Point>>>& gestures, int e) {
		//wxlogmessage("%d", gestures.size());
	    
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine randend(seed);
		vector<splitDataStruct> totaldata;
		srand(time(nullptr));
		bool flag = false;
	
		for (auto& it : gestures) {
			splitDataStruct temp;
			temp.name = it.first;
			temp.gesturePoints = it.second;
			totaldata.push_back(temp);
			temp = splitDataStruct();
		}
		/*for (auto& it : totaldata) {
			wxlogmessage("gestures name %s ", it.name);
			for (auto& it1 : it.gesturepoints) {
				for (auto& it2 : it1) {
					wxlogmessage("%f%f",it2.x,it2.y);
				}
			}
		}*/
		//i + 16 | i = 0; i < 15;
		auto& pointsvec = totaldata;
		//wxLogMessage("size: %d", totaldata.size());
		vector<int> indices(pointsvec.size());
		iota(indices.begin(), indices.end(), 0);
		shuffle(indices.begin(), indices.end(), randend);
		vector<int> trainingIndices(10);
		iota(trainingIndices.begin(), trainingIndices.end(), 0);
		shuffle(trainingIndices.begin(), trainingIndices.end(), randend);
		/*wxlogmessage("map size %d", pointsvec.size());
		for (auto& it : indices) {
			wxlogmessage("indices in the vector %d", it);
		}*/

			splitDataStruct testdata;
			vector<splitDataStruct> trainingdata;
			// check that the size of the gesture vector is at least 
			int index = indices[0];
			testdata.name = pointsvec[index].name;
			testdata.gesturePoints = pointsvec[index].gesturePoints;
			for (int j = 1; j <= totaldata.size(); j+=10) {
				splitDataStruct temp;
				for (int i = 1; i <= e; i++) {
					int index = trainingIndices[i];
					if (j + index != 160) {
						temp.name = totaldata[j + index].name;
						temp.gesturePoints = totaldata[j + index].gesturePoints;
						shuffle(trainingIndices.begin(), trainingIndices.end(), randend);
						trainingdata.push_back(temp);
					}
					else {
						temp.name = totaldata[159].name;
						temp.gesturePoints = totaldata[159].gesturePoints;
						shuffle(trainingIndices.begin(), trainingIndices.end(), randend);
						trainingdata.push_back(temp);
					}
				}
				temp = splitDataStruct();
			}
			
			//for (auto& it : totaldata) {
			//	wxLogMessage("gestures name %s ", it.name);
			//	/*for (auto& it1 : it.gesturepoints) {
			//		for (auto& it2 : it1) {
			//			wxlogmessage("%f%f", it2.x, it2.y);
			//		}
			//	}*/
			//}

			
		
		return make_pair(testdata, trainingdata);
	}

};
struct CollectDataStruct {
	string filename;
	string subject;
	string date;
	string TimeOfDay;
	string speed;

};

class CollectData {
public:
	void savePointsToXml(const vector<vector<Point>>& points, const std::string& filename, string& gesture)
	{
		// Create the XML document

		tinyxml2::XMLDocument doc;
		tinyxml2::XMLDeclaration* decl = doc.NewDeclaration();
		doc.InsertFirstChild(decl);

		// Create the root element
		tinyxml2::XMLElement* root = doc.NewElement("Gesture");
		root->SetAttribute("Name", gesture.c_str());
		root->SetAttribute("Subject", "2");
		root->SetAttribute("InputType", "Touch");
		/*root->SetAttribute("Speed", "medium");
		root->SetAttribute("Number", "1");*/
		root->SetAttribute("NumPts", points.size()*points[0].size());
		
	/*	root->SetAttribute("Millseconds", "1268");
		root->SetAttribute("AppName", "Gestures");
		root->SetAttribute("AppVer", "3.5.0.0");
		root->SetAttribute("Date", "Monday, March 05, 2007");*/
		root->SetAttribute("data-type", "simple");
		doc.InsertEndChild(root);
		// Add each point as a child element
		//wxLogMessage("storing\n");
		int strokeIndex = 1;
		while (strokeIndex <= points.size()) {
			tinyxml2::XMLElement* strokeElement = doc.NewElement("Stroke");
			strokeElement->SetAttribute("index", strokeIndex);
			root->InsertEndChild(strokeElement);
			for (const auto& strokes : points)
			{	
				for (const auto& p : strokes) {
					tinyxml2::XMLElement* pointElement = doc.NewElement("Point");
					pointElement->SetAttribute("X", p.x);
					pointElement->SetAttribute("Y", p.y);
					pointElement->SetAttribute("T", p.t);
					pointElement->SetAttribute("Pressure", "128");
					strokeElement->InsertEndChild(pointElement);
				}
			}
			
			strokeIndex++;
		}
		// Create the folder if it doesn't exist
		std::filesystem::path folderPath("user_inputs");
		if (!std::filesystem::is_directory(folderPath))
		{
			if (!std::filesystem::create_directory(folderPath))
			{
				std::cerr << "Failed to create directory: " << folderPath << std::endl;
				return;
			}
		}
		//wxLogMessage("writing\n");
		// Save the document to file in the folder
		std::filesystem::path filePath = folderPath / filename;
		if (doc.SaveFile(filePath.string().c_str()) != tinyxml2::XML_SUCCESS)
		{
			std::cerr << "Failed to save XML file: " << filePath << std::endl;
		}
		// Save the document to file
	}

};







class MyCanvas : public wxWindow
{
public:
	//Canvas for drawing 
	MyCanvas(wxWindow* parent) : wxWindow(parent, wxID_ANY)
	{

		m_prompt = new wxStaticText(this, wxID_ANY, "click submit to begin");
		//m_prompt = new wxStaticText(this, wxID_ANY, "Draw a {Gesture}");
		m_output = new wxStaticText(this, wxID_ANY, "");
		m_counter = new wxStaticText(this, wxID_ANY, "0 / 160");
		//total_counter = new wxStaticText(this, wxID_ANY, "Total Counter: 0/160");

		m_clearButton = new wxButton(this, wxID_ANY, "Clear");
		m_submitButton = new wxButton(this, wxID_ANY, "Submit");

		wxFont font(10, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
		m_font = font;
		//m_output->SetFont(m_font);

		//total_counter->SetFont(m_font);
		m_prompt->SetFont(m_font);
		m_output->SetFont(m_font);
		m_counter->SetFont(m_font);

		m_submitButton->SetBackgroundColour(wxColour(0, 255, 0));
		m_clearButton->SetBackgroundColour(wxColour(255, 0, 0));

		//Main Sizer 
		wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);

		//sizer->Add(m_output, 0, wxALIGN_CENTER_HORIZONTAL | wxTOP, 0);
		//Text Sizer
		wxBoxSizer* textSizer = new wxBoxSizer(wxHORIZONTAL);
		textSizer->Add(m_prompt, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 0);
		textSizer->Add(m_output, 0, wxALIGN_CENTER_VERTICAL | wxTOP, 0);
		//textSizer->Add(total_counter, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 0); // add counter to sizer

		//Button Sizer
		wxBoxSizer* buttonsSizer = new wxBoxSizer(wxHORIZONTAL);
		buttonsSizer->Add(m_clearButton, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 20);
		buttonsSizer->Add(m_submitButton, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 20);
		//creating space from text
		sizer->Add(textSizer, 0, wxALIGN_CENTRE_HORIZONTAL | wxTop, 10);
		sizer->AddStretchSpacer(1);
		sizer->Add(buttonsSizer, 0, wxALIGN_CENTRE_HORIZONTAL | wxBottom, 10);



		//sizer->Add(m_clearButton, 0, wxALIGN_CENTER_HORIZONTAL | wxBottom, 0);
		//sizer->Add(m_submitButton, 0, wxALIGN_CENTER_HORIZONTAL | wxBottom, 0);
		SetSizerAndFit(sizer);

		Connect(wxEVT_PAINT, wxPaintEventHandler(MyCanvas::OnPaint));
		Connect(wxEVT_LEFT_DOWN, wxMouseEventHandler(MyCanvas::OnLeftDown));
		Connect(wxEVT_LEFT_UP, wxMouseEventHandler(MyCanvas::OnLeftUp));
		Connect(wxEVT_MOTION, wxMouseEventHandler(MyCanvas::OnMotion));
		Connect(m_clearButton->GetId(), wxEVT_BUTTON, wxCommandEventHandler(MyCanvas::OnClearButtonClick));
		Connect(m_submitButton->GetId(), wxEVT_BUTTON, wxCommandEventHandler(MyCanvas::OnSubmitButtonClick));
		OfflineRecognizer a;
		a.recognizeOfflineData();
		//a.preProcessOfflineData();
	}
private:
	//event handler is called when the canvas is repainted
	wxStaticText* m_prompt;
	wxStaticText* m_output;
	wxStaticText* m_counter;
	wxStaticText* total_counter;
	wxButton* m_clearButton;
	bool m_isDrawing = false;
	wxButton* m_submitButton;
	vector<wxPoint> m_points;
	vector<Point> tpoints;
	vector<vector<Point>> strokes;
	wxFont m_font;
	vector<Point> points;
	unordered_map<string, int> mp;
	int counter = -1;
	int current_gesture = 0;
	
	vector<string> labelList = { "arrowhead", "asterisk", "D", "exclamation_point", "five_point_star", "H", "half_note", "I","line", "N", "null", "P", "pitchfork","six_point_star","T", "X" };
	//vector<string> labelList = { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P" };

	
	CollectData dataCollection;
	
	void OnPaint(wxPaintEvent& event)
	{
		wxPaintDC dc(this);
		wxGraphicsContext* gc = wxGraphicsContext::Create(dc);

		static vector<wxColour> colors = { wxColour(255, 0, 0), wxColour(0, 255, 0), wxColour(0, 0, 255) };
		static int colorIndex = 0;

		gc->SetPen(wxPen(colors[colorIndex], 5));

		size_t i = 1;
		while (i < m_points.size()) {
			while (i < m_points.size() && m_points[i - 1] != wxPoint(-1, -1) && m_points[i] != wxPoint(-1, -1)) {
				dc.DrawLine(m_points[i - 1].x, m_points[i - 1].y, m_points[i].x, m_points[i].y);
				i++;
			}
			i++;
		}

		if (!m_isDrawing) {
			colorIndex = (colorIndex + 1) % colors.size();
			delete gc;
		}
	}

		void OnSubmitButtonClick(wxCommandEvent & event) {
			if (counter == -1) {
				//current_gesture++;
				counter = 0;
				wxString next_gesture = wxString::Format("Draw a %s", labelList[current_gesture]);
				m_prompt->SetLabel(next_gesture);
				return;
			}
			if (counter < 160 && mp[labelList[current_gesture]] != 10) {
				vector<vector<Point>> strokes;
				std::vector<Point> currentStroke;
				
				for (size_t i = 1; i < tpoints.size(); i++) {
					Point t = tpoints[i];
					//const wxPoint& p = m_points[i];

					// Check if this is the end of a stroke
					//wxLogMessage("POINT P: (%f,%f)", tpoints[i].x, tpoints[i].y);
					if (t.x == -1.0 && t.y == -1.0 || i == tpoints.size() - 1) {
						//wxLogMessage(" Inside POINT P: (%f,%f)", t.x, t.y);
						// Add the last point if it hasn't already been added
						//if (i != m_points.size() - 1) {
					
						if (t.x != -1 && t.y != -1){
							currentStroke.emplace_back(t.x, t.y,t.t);
						}

						// Add the current stroke to the list of strokes and start a new stroke
						if (currentStroke.size() > 2)
							strokes.push_back(currentStroke);
						currentStroke.clear();
					}
					else {
						currentStroke.emplace_back(t.x, t.y, t.t);
					}
				}
				
				mp[labelList[current_gesture]]++;
				string name = labelList[current_gesture] + to_string(mp[labelList[current_gesture]]) + ".xml";
				dataCollection.savePointsToXml(strokes, name,labelList[current_gesture]);
				counter++;
				/*current_gesture++ % labelList.size()-1;*/
				wxString next_gesture = wxString::Format("Draw a %s", labelList[current_gesture]);
				m_prompt->SetLabel(next_gesture);

				//[VIKRANTH]: gesture is drawn so you can now save it in the file with the name string name = labelList[current_gesture] + to_string(mp[labelList[current_gesture]])
				//and the vector is already ready so you just need to push it. 
				//data structure map<string(gesture name), map<string(useName), vector<Point>(storing the points)>>>
			}
			else if (counter < 160 && mp[labelList[current_gesture]] == 10) {
				//current_gesture++;
				/*for (auto label : labelList) {
					wxLogMessage("Before %s\n", label);
				}*/
				/*auto it = labelList.begin() + current_gesture;
				labelList.erase(it);*/
				current_gesture++;
				wxString next_gesture = wxString::Format("Draw a %s", labelList[current_gesture]);
				m_prompt->SetLabel(next_gesture);
				/*for (auto label : labelList) {
					wxLogMessage("after %s\n", label);
				}*/
				/*wxString next_gesture = wxString::Format("Draw a %s", labelList[current_gesture]);
				m_prompt->SetLabel(next_gesture);*/
				//continue;
			}
			else {
				//TO DO VIKRANTH : save these points into an xml file, I think I have also laid out the flow. 
				m_prompt->SetLabel("All gestures complete");
				for (auto label : mp) {
					wxLogMessage("Gesture %s: %d", label.first, label.second);
				}
			}



			//vector<vector<Point>> strokes;
			//std::vector<Point> currentStroke;
			//for (size_t i = 1; i < m_points.size(); i++) {
			//	const wxPoint& p = m_points[i];

			//	// Check if this is the end of a stroke
			//	//wxLogMessage("POINT P: (%d,%d)", m_points[i].x,m_points[i].y);
			//	if (p.x== -1 && p.y== -1 || i == m_points.size() - 1) {
			//		//wxLogMessage(" Inside POINT P: (%d,%d)", p.x, p.y);
			//		// Add the last point if it hasn't already been added
			//		//if (i != m_points.size() - 1) {
			//			if(p.x!=-1 &&p.y!=-1)
			//			currentStroke.emplace_back(p.x, p.y);
			//		//}

			//		// Add the current stroke to the list of strokes and start a new stroke
			//		if(currentStroke.size()>2)
			//		strokes.push_back(currentStroke);
			//		currentStroke.clear();
			//	}
			//	else {
			//		currentStroke.emplace_back(p.x, p.y);
			//	}
			//}

			//MGestureRecognizer GR(false); // Instantiate your recognizer object
			//Result res = GR.Recognize(strokes, false, false, false);
			//wxString outputStr = wxString::Format("Result: %s (%f) in %d ms", res.Name, res.Score, res.Time);
		/*	m_prompt->SetLabel(outputStr);*/
			tpoints.clear();
			m_points.clear();
			Refresh();
			wxString counterStr = wxString::Format("Counter: %d", counter);
			m_counter->SetLabel(counterStr + "/160");


		}
		void OnClearButtonClick(wxCommandEvent& event)
		{
			m_points.clear();
			Refresh();
		}
	//event handler is called when the left mouse button is pressed and pushes position to m_points
	void OnLeftDown(wxMouseEvent& event)
	{
		if (!m_isDrawing) {
			m_isDrawing = true;
			m_points.push_back(wxPoint(-1, -1));
			tpoints.push_back(Point(-1, -1));
		}
		//wxDateTime now = wxDateTime::Now();
		auto now = std::chrono::system_clock::now();
		auto duration = now.time_since_epoch();
		auto seconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
		auto time_since_epoch = static_cast<std::time_t>(seconds.count());
		//long long milliseconds = now.GetTicks();
		//int seconds = (int)(milliseconds / 1000);
		wxPoint p = event.GetPosition();
		tpoints.push_back(Point(p.x, p.y, abs(time_since_epoch)));
		tpoints.push_back(Point(-1, -1));
		m_points.push_back(event.GetPosition());
		CaptureMouse();
	}

	//event handler is called when the left mouse button is released

	void OnLeftUp(wxMouseEvent& event)
	{
		auto now = std::chrono::system_clock::now();
		auto duration = now.time_since_epoch();
		auto seconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
		auto time_since_epoch = static_cast<std::time_t>(seconds.count());
		wxPoint p = event.GetPosition();
		tpoints.push_back(Point(p.x, p.y, abs(time_since_epoch)));
		tpoints.push_back(Point(-1, -1));
		m_points.push_back(event.GetPosition());
		ReleaseMouse();
		m_points.push_back(wxPoint(-1, -1));
		m_isDrawing = false;
	}
	//void OnLeftUp(wxMouseEvent& event)
	//{
	//	if (m_isDrawing)
	//	{
	//		
	//		m_points.push_back(event.GetPosition());
	//		m_points.push_back(wxPoint(-1, -1));
	//		m_isDrawing = false;
	//		ReleaseMouse();
	//	}
	//}

	//function event handler is called when the mouse is moved
	void OnMotion(wxMouseEvent& event)
	{
		if (event.Dragging() && event.LeftIsDown() && m_isDrawing)
		{
			auto now = std::chrono::system_clock::now();
			auto duration = now.time_since_epoch();
			auto seconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
			auto time_since_epoch = static_cast<std::time_t>(seconds.count());
			wxPoint p = event.GetPosition();
			tpoints.push_back(Point(p.x, p.y, abs(time_since_epoch)));
			//m_output->SetLabel("Recording unistroke...");
			m_points.push_back(event.GetPosition());
			Refresh();
		}
	}
	

};

class MyApp : public wxApp
{
public:
	virtual bool OnInit()
	{
		wxFrame* frame = new wxFrame(NULL, wxID_ANY, "$1 Recognizer");
		MyCanvas* canvas = new MyCanvas(frame);
		frame->Show();


		return true;
	}
};


//creates an instance of the MyApp
wxIMPLEMENT_APP(MyApp);