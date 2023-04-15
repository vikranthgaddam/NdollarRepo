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
using namespace std;

class Point {//Defining point struct 
public:
	double x;
	double y;
	Point(double x, double y) {
		this->x = x;
		this->y = y;
	}
	Point() {}
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




class MyCanvas : public wxWindow
{
public:
	//Canvas for drawing 
	MyCanvas(wxWindow* parent) : wxWindow(parent, wxID_ANY)
	{

		m_output = new wxStaticText(this, wxID_ANY, "");
		m_submitButton = new wxButton(this, wxID_ANY, "Submit");
		wxFont font(12, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
		m_font = font;
		m_output->SetFont(m_font);

		m_submitButton->SetBackgroundColour(wxColour(0, 255, 0));
		wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
		//Text Sizer
		wxBoxSizer* textSizer = new wxBoxSizer(wxHORIZONTAL);
		textSizer->Add(m_output, 0, wxALIGN_CENTER_VERTICAL | wxLeft, 0);
		sizer->Add(textSizer, 0, wxALIGN_CENTRE_HORIZONTAL | wxTop, 10);
		sizer->AddStretchSpacer(1);

		wxBoxSizer* buttonsSizer = new wxBoxSizer(wxHORIZONTAL);
		buttonsSizer->Add(m_submitButton, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 20);

		sizer->Add(buttonsSizer, 0, wxALIGN_CENTRE_HORIZONTAL | wxBottom, 10);

		SetSizerAndFit(sizer);

		Connect(wxEVT_PAINT, wxPaintEventHandler(MyCanvas::OnPaint));
		Connect(wxEVT_LEFT_DOWN, wxMouseEventHandler(MyCanvas::OnLeftDown));
		Connect(wxEVT_LEFT_UP, wxMouseEventHandler(MyCanvas::OnLeftUp));
		Connect(wxEVT_MOTION, wxMouseEventHandler(MyCanvas::OnMotion));
		Connect(m_submitButton->GetId(), wxEVT_BUTTON, wxCommandEventHandler(MyCanvas::OnSubmitButtonClick));
	}
private:
	//event handler is called when the canvas is repainted
	bool m_isDrawing = false;
	wxButton* m_submitButton;
	vector<wxPoint> m_points;
	vector<vector<wxPoint>> listOf;
	vector<vector<Point>> strokes;
	wxFont m_font;
	vector<Point> points;
	wxStaticText* m_output;
	void OnPaint(wxPaintEvent& event)
	{
		wxPaintDC dc(this);
		wxGraphicsContext* gc = wxGraphicsContext::Create(dc);

		gc->SetPen(wxPen(*wxRED, 5));

		size_t i = 1;
		while (i < m_points.size()) {
			while (i < m_points.size() && m_points[i - 1] != wxPoint(-1, -1) && m_points[i] != wxPoint(-1, -1)) {
				dc.DrawLine(m_points[i - 1].x, m_points[i - 1].y, m_points[i].x, m_points[i].y);
				i++;
			}
			i++;
		}

		delete gc;
	}

		void OnSubmitButtonClick(wxCommandEvent & event) {

			vector<vector<Point>> strokes;
			std::vector<Point> currentStroke;
			for (size_t i = 1; i < m_points.size(); i++) {
				const wxPoint& p = m_points[i];

				// Check if this is the end of a stroke
				//wxLogMessage("POINT P: (%d,%d)", m_points[i].x,m_points[i].y);
				if (p.x== -1 && p.y== -1 || i == m_points.size() - 1) {
					//wxLogMessage(" Inside POINT P: (%d,%d)", p.x, p.y);
					// Add the last point if it hasn't already been added
					//if (i != m_points.size() - 1) {
						if(p.x!=-1 &&p.y!=-1)
						currentStroke.emplace_back(p.x, p.y);
					//}

					// Add the current stroke to the list of strokes and start a new stroke
					if(currentStroke.size()>2)
					strokes.push_back(currentStroke);
					currentStroke.clear();
				}
				else {
					currentStroke.emplace_back(p.x, p.y);
				}
			}
			wxLogMessage("Strokes:");
			for (size_t i = 0; i < strokes.size(); i++) {
				wxLogMessage("Stroke %d:", i);
				const vector<Point>& points = strokes[i];
				for (size_t j = 0; j < points.size(); j++) {
					wxLogMessage("    (%f, %f)", points[j].x, points[j].y);
				}
			}
			MGestureRecognizer GR(false); // Instantiate your recognizer object
			Result res = GR.Recognize(strokes, false, false, false);
			wxString outputStr = wxString::Format("Result: %s (%f) in %d ms", res.Name, res.Score, res.Time);
			m_output->SetLabel(outputStr);
			m_points.clear();
			Refresh();
		}

	//event handler is called when the left mouse button is pressed and pushes position to m_points
	void OnLeftDown(wxMouseEvent& event)
	{
		if (!m_isDrawing) {
			m_isDrawing = true;
			m_points.push_back(wxPoint(-1, -1));
		}
		m_points.push_back(event.GetPosition());
		CaptureMouse();
	}

	//event handler is called when the left mouse button is released
	void OnLeftUp(wxMouseEvent& event)
	{
		m_points.push_back(event.GetPosition());
		ReleaseMouse();
		m_points.push_back(wxPoint(-1, -1));
		m_isDrawing = false;
	}

	//function event handler is called when the mouse is moved
	void OnMotion(wxMouseEvent& event)
	{
		if (event.Dragging() && event.LeftIsDown() && m_isDrawing)
		{
			m_output->SetLabel("Recording unistroke...");
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