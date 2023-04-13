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

using namespace std;

class Point {//Defining point struct 
public:
	double x;
	double y;
	Point(double x, double y) {
		this->x = x;
		this->y = y;
	}
};

struct Rect { //Defining rectangle for bounding box
	double X;
	double Y;
	double Width;
	double Height;

	Rect(double x, double y, double width, double height)
		: X(x), Y(y), Width(width), Height(height) {}
};

const int NumUnistrokes = 16;//Tempates loaded
const int NumPoints = 64;
const double SquareSize = 250.0;
const Point Origin = { 0, 0 };
const double Diagonal = sqrt(SquareSize * SquareSize + SquareSize * SquareSize);//Diagonal length
const double HalfDiagonal = 0.5 * Diagonal;
const double AngleRange = 45.0 * M_PI / 180.0;
const double AnglePrecision = 2.0 * M_PI / 180.0;
const double Phi = 0.5 * (-1.0 + sqrt(5.0));


struct Unistroke {//This is where preprocessing happens for points
public:
	std::string name;
	std::vector<Point> points;
	std::vector<double> vectorizedPoints;


	Unistroke(std::string& name, std::vector<Point>& points) //Constructor where the point are preprocessed and stored
	{
		//Note 1.Resample 2.Rotate By needs indicative angle 3.Scale 4.Transalate 5.For Protractor -- Vectorize
		this->name = name;
		this->points = Resample(points, NumPoints);//1
		double radians = IndicativeAngle(this->points);//2
		this->points = RotateBy(this->points, -radians);//3
		this->points = ScaleTo(this->points, SquareSize);//4
		this->points = TranslateTo(this->points, Origin);//5
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

	std::vector<Point> ScaleTo(std::vector<Point> points, double size)
	{
		Rect B = BoundingBox(points);
		std::vector<Point> newpoints;
		for (int i = 0; i < points.size(); i++)
		{
			double qx = points[i].x * (size / B.Width);
			double qy = points[i].y * (size / B.Height);
			newpoints.push_back(Point(qx, qy));
		}//TODO Raghav
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
class GestureRecognizer {

public:
	GestureRecognizer() {
		//Intialize a empty vector with type Unistroke
		this->Unistrokes = std::vector<Unistroke>();
		//Note:Need to add points here 
		std::vector<Point> trianglePoints = {
		Point(137, 139), Point(135, 141), Point(133, 144), Point(132, 146),
		Point(130, 149), Point(128, 151), Point(126, 155), Point(123, 160),
		Point(120, 166), Point(116, 171), Point(112, 177), Point(107, 183),
		Point(102, 188), Point(100, 191), Point(95, 195), Point(90, 199),
		Point(86, 203), Point(82, 206), Point(80, 209), Point(75, 213),
		Point(73, 213), Point(70, 216), Point(67, 219), Point(64, 221),
		Point(61, 223), Point(60, 225), Point(62, 226), Point(65, 225),
		Point(67, 226), Point(74, 226), Point(77, 227), Point(85, 229),
		Point(91, 230), Point(99, 231), Point(108, 232), Point(116, 233),
		Point(125, 233), Point(134, 234), Point(145, 233), Point(153, 232),
		Point(160, 233), Point(170, 234), Point(177, 235), Point(179, 236),
		Point(186, 237), Point(193, 238), Point(198, 239), Point(200, 237),
		Point(202, 239), Point(204, 238), Point(206, 234), Point(205, 230),
		Point(202, 222), Point(197, 216), Point(192, 207), Point(186, 198),
		Point(179, 189), Point(174, 183), Point(170, 178), Point(164, 171),
		Point(161, 168), Point(154, 160), Point(148, 155), Point(143, 150),
		Point(138, 148), Point(136, 148)
		};
		std::vector<Point> xPoints = { Point(87, 142),  Point(89, 145),  Point(91, 148),  Point(93, 151)
			,  Point(96, 155),  Point(98, 157),  Point(100, 160),  Point(102, 162),  Point(106, 167),
			 Point(108, 169),  Point(110, 171),  Point(115, 177),  Point(119, 183), Point(123, 189),
			 Point(127, 193),  Point(129, 196),  Point(133, 200),  Point(137, 206),  Point(140, 209),
			 Point(143, 212),  Point(146, 215),  Point(151, 220),  Point(153, 222),  Point(155, 223),
			 Point(157, 225),  Point(158, 223),  Point(157, 218),  Point(155, 211),  Point(154, 208),
			 Point(152, 200),  Point(150, 189),  Point(148, 179),  Point(147, 170),  Point(147, 158),
			 Point(147, 148),  Point(147, 141),  Point(147, 136),  Point(144, 135),  Point(142, 137),
			 Point(140, 139),  Point(135, 145),  Point(131, 152),  Point(124, 163),  Point(116, 177),
			 Point(108, 191),  Point(100, 206),  Point(94, 217),  Point(91, 222),  Point(89, 225),  Point(87, 226),  Point(87, 224) };
		std::vector<Point> xrectanglePoints = { Point(78, 149),  Point(78, 153),  Point(78, 157),  Point(78, 160),  Point(79, 162),  Point(79, 164),  Point(79, 167),
			Point(79, 169), Point(79, 173),  Point(79, 178),  Point(79, 183),  Point(80, 189),  Point(80, 193),  Point(80, 198),  Point(80, 202),
			 Point(81, 208),  Point(81, 210),  Point(81, 216),  Point(82, 222),  Point(82, 224),  Point(82, 227),  Point(83, 229),  Point(83, 231),
			 Point(85, 230),  Point(88, 232),  Point(90, 233),  Point(92, 232),  Point(94, 233),  Point(99, 232),  Point(102, 233),  Point(106, 233),
			 Point(109, 234),  Point(117, 235),  Point(123, 236),  Point(126, 236),  Point(135, 237),  Point(142, 238),  Point(145, 238),  Point(152, 238),
			 Point(154, 239),  Point(165, 238),  Point(174, 237),  Point(179, 236),  Point(186, 235),  Point(191, 235),  Point(195, 233),  Point(197, 233),
			 Point(200, 233),  Point(201, 235),  Point(201, 233),  Point(199, 231),  Point(198, 226),  Point(198, 220),  Point(196, 207),  Point(195, 195),
			 Point(195, 181),  Point(195, 173),  Point(195, 163),  Point(194, 155),  Point(192, 145),  Point(192, 143),  Point(192, 138),  Point(191, 135),
			 Point(191, 133),  Point(191, 130),  Point(190, 128),  Point(188, 129),  Point(186, 129),  Point(181, 132),  Point(173, 131),  Point(162, 131),
			 Point(151, 132),  Point(149, 132),  Point(138, 132),  Point(136, 132),  Point(122, 131),  Point(120, 131),  Point(109, 130),  Point(107, 130),
			 Point(90, 132),  Point(81, 133),  Point(76, 133) };
		std::vector<Point> circle = { Point(127, 141),  Point(124, 140),  Point(120, 139),  Point(118, 139),  Point(116, 139),
			 Point(111, 140),  Point(109, 141),  Point(104, 144),  Point(100, 147),  Point(96, 152),  Point(93, 157),  Point(90, 163),
			 Point(87, 169),  Point(85, 175),  Point(83, 181),  Point(82, 190),  Point(82, 195),  Point(83, 200),  Point(84, 205),
			 Point(88, 213),  Point(91, 216),  Point(96, 219),  Point(103, 222),  Point(108, 224),  Point(111, 224),  Point(120, 224),
			 Point(133, 223),  Point(142, 222),  Point(152, 218),  Point(160, 214),  Point(167, 210),  Point(173, 204),  Point(178, 198),
			 Point(179, 196),  Point(182, 188),  Point(182, 177),  Point(178, 167),  Point(170, 150),  Point(163, 138),  Point(152, 130),
			 Point(143, 129),  Point(140, 131),  Point(129, 136),  Point(126, 139) };
		std::vector<Point> checkPoints = { Point(91, 185), Point(93, 185), Point(95, 185), Point(97, 185), Point(100, 188), Point(102, 189), Point(104, 190), Point(106, 193), Point(108, 195), Point(110, 198), Point(112, 201), Point(114, 204), Point(115, 207), Point(117, 210), Point(118, 212), Point(120, 214), Point(121, 217), Point(122, 219), Point(123, 222), Point(124, 224), Point(126, 226), Point(127, 229), Point(129, 231), Point(130, 233), Point(129, 231), Point(129, 228), Point(129, 226), Point(129, 224), Point(129, 221), Point(129, 218), Point(129, 212), Point(129, 208), Point(130, 198), Point(132, 189), Point(134, 182), Point(137, 173), Point(143, 164), Point(147, 157), Point(151, 151), Point(155, 144), Point(161, 137), Point(165, 131), Point(171, 122), Point(174, 118), Point(176, 114), Point(177, 112), Point(177, 114), Point(175, 116), Point(173, 118) };
		std::vector<Point> caret = { Point(79, 245), Point(79, 242), Point(79, 239), Point(80, 237), Point(80, 234), Point(81, 232), Point(82, 230), Point(84, 224),
			Point(86, 220), Point(86, 218), Point(87, 216), Point(88, 213), Point(90, 207), Point(91, 202), Point(92, 200), Point(93, 194), Point(94, 192), Point(96, 189),
			Point(97, 186), Point(100, 179), Point(102, 173), Point(105, 165), Point(107, 160), Point(109, 158), Point(112, 151), Point(115, 144), Point(117, 139),
			Point(119, 136), Point(119, 134), Point(120, 132), Point(121, 129), Point(122, 127), Point(124, 125), Point(126, 124), Point(129, 125), Point(131, 127),
			Point(132, 130), Point(136, 139), Point(141, 154), Point(145, 166), Point(151, 182), Point(156, 193), Point(157, 196), Point(161, 209), Point(162, 211),
			Point(167, 223), Point(169, 229), Point(170, 231), Point(173, 237), Point(176, 242), Point(177, 244), Point(179, 250), Point(181, 255), Point(182, 257) };
		std::vector<Point>zigZag = { Point(307, 216), Point(333, 186), Point(356, 215), Point(375, 186), Point(399, 216), Point(418, 186) };
		std::vector<Point>arrow = { Point(68, 222), Point(70, 220), Point(73, 218), Point(75, 217), Point(77, 215), Point(80, 213), Point(82, 212), Point(84, 210), Point(87, 209), Point(89, 208), Point(92, 206), Point(95, 204), Point(101, 201), Point(106, 198), Point(112, 194), Point(118, 191), Point(124, 187), Point(127, 186), Point(132, 183), Point(138, 181), Point(141, 180), Point(146, 178), Point(154, 173), Point(159, 171), Point(161, 170), Point(166, 167), Point(168, 167), Point(171, 166), Point(174, 164), Point(177, 162), Point(180, 160), Point(182, 158), Point(183, 156), Point(181, 154), Point(178, 153), Point(171, 153), Point(164, 153), Point(160, 153), Point(150, 154), Point(147, 155), Point(141, 157), Point(137, 158), Point(135, 158), Point(137, 158), Point(140, 157), Point(143, 156), Point(151, 154), Point(160, 152), Point(170, 149), Point(179, 147), Point(185, 145), Point(192, 144), Point(196, 144), Point(198, 144), Point(200, 144), Point(201, 147), Point(199, 149), Point(194, 157), Point(191, 160), Point(186, 167), Point(180, 176), Point(177, 179), Point(171, 187), Point(169, 189), Point(165, 194), Point(164, 196)
		};

		std::vector<Point> leftsquarebracket = { Point(140, 124), Point(138, 123), Point(135, 122), Point(133, 123), Point(130, 123), Point(128, 124), Point(125, 125), Point(122, 124), Point(120, 124), Point(118, 124), Point(116, 125), Point(113, 125), Point(111, 125), Point(108, 124), Point(106, 125), Point(104, 125), Point(102, 124), Point(100, 123), Point(98, 123), Point(95, 124), Point(93, 123), Point(90, 124), Point(88, 124), Point(85, 125), Point(83, 126), Point(81, 127), Point(81, 129), Point(82, 131), Point(82, 134), Point(83, 138), Point(84, 141), Point(84, 144), Point(85, 148), Point(85, 151), Point(86, 156), Point(86, 160), Point(86, 164), Point(86, 168), Point(87, 171), Point(87, 175), Point(87, 179), Point(87, 182), Point(87, 186), Point(88, 188), Point(88, 195), Point(88, 198), Point(88, 201), Point(88, 207), Point(89, 211), Point(89, 213), Point(89, 217), Point(89, 222), Point(88, 225), Point(88, 229), Point(88, 231), Point(88, 233), Point(88, 235), Point(89, 237), Point(89, 240), Point(89, 242), Point(91, 241), Point(94, 241), Point(96, 240), Point(98, 239), Point(105, 240), Point(109, 240), Point(113, 239), Point(116, 240), Point(121, 239), Point(130, 240), Point(136, 237), Point(139, 237), Point(144, 238), Point(151, 237), Point(157, 236), Point(159, 237) };
		std::vector<Point> rightsquarebracket = { Point(112, 138), Point(112, 136), Point(115, 136), Point(118, 137), Point(120, 136), Point(123, 136), Point(125, 136), Point(128, 136), Point(131, 136), Point(134, 135), Point(137, 135), Point(140, 134), Point(143, 133), Point(145, 132), Point(147, 132), Point(149, 132), Point(152, 132), Point(153, 134), Point(154, 137), Point(155, 141), Point(156, 144), Point(157, 152), Point(158, 161), Point(160, 170), Point(162, 182), Point(164, 192), Point(166, 200), Point(167, 209), Point(168, 214), Point(168, 216), Point(169, 221), Point(169, 223), Point(169, 228), Point(169, 231), Point(166, 233), Point(164, 234), Point(161, 235), Point(155, 236), Point(147, 235), Point(140, 233), Point(131, 233), Point(124, 233), Point(117, 235), Point(114, 238), Point(112, 238) };
		std::vector<Point>v = { Point(89, 164), Point(90, 162), Point(92, 162), Point(94, 164), Point(95, 166), Point(96, 169), Point(97, 171), Point(99, 175), Point(101, 178), Point(103, 182), Point(106, 189), Point(108, 194), Point(111, 199), Point(114, 204), Point(117, 209), Point(119, 214), Point(122, 218), Point(124, 222), Point(126, 225), Point(128, 228), Point(130, 229), Point(133, 233), Point(134, 236), Point(136, 239), Point(138, 240), Point(139, 242), Point(140, 244), Point(142, 242), Point(142, 240), Point(142, 237), Point(143, 235), Point(143, 233), Point(145, 229), Point(146, 226), Point(148, 217), Point(149, 208), Point(149, 205), Point(151, 196), Point(151, 193), Point(153, 182), Point(155, 172), Point(157, 165), Point(159, 160), Point(162, 155), Point(164, 150), Point(165, 148), Point(166, 146) };
		std::vector<Point>deletePoints = { Point(123, 129), Point(123, 131), Point(124, 133), Point(125, 136), Point(127, 140), Point(129, 142), Point(133, 148), Point(137, 154), Point(143, 158), Point(145, 161), Point(148, 164), Point(153, 170), Point(158, 176), Point(160, 178), Point(164, 183), Point(168, 188), Point(171, 191), Point(175, 196), Point(178, 200), Point(180, 202), Point(181, 205), Point(184, 208), Point(186, 210), Point(187, 213), Point(188, 215), Point(186, 212), Point(183, 211), Point(177, 208), Point(169, 206), Point(162, 205), Point(154, 207), Point(145, 209), Point(137, 210), Point(129, 214), Point(122, 217), Point(118, 218), Point(111, 221), Point(109, 222), Point(110, 219), Point(112, 217), Point(118, 209), Point(120, 207), Point(128, 196), Point(135, 187), Point(138, 183), Point(148, 167), Point(157, 153), Point(163, 145), Point(165, 142), Point(172, 133), Point(177, 127), Point(179, 127), Point(180, 125) };
		std::vector<Point> leftCB = { Point(150, 116), Point(147, 117), Point(145, 116), Point(142, 116), Point(139, 117), Point(136, 117), Point(133, 118), Point(129, 121), Point(126, 122), Point(123, 123), Point(120, 125), Point(118, 127), Point(115, 128), Point(113, 129), Point(112, 131), Point(113, 134), Point(115, 134), Point(117, 135), Point(120, 135), Point(123, 137), Point(126, 138), Point(129, 140), Point(135, 143), Point(137, 144), Point(139, 147), Point(141, 149), Point(140, 152), Point(139, 155), Point(134, 159), Point(131, 161), Point(124, 166), Point(121, 166), Point(117, 166), Point(114, 167), Point(112, 166), Point(114, 164), Point(116, 163), Point(118, 163), Point(120, 162), Point(122, 163), Point(125, 164), Point(127, 165), Point(129, 166), Point(130, 168), Point(129, 171), Point(127, 175), Point(125, 179), Point(123, 184), Point(121, 190), Point(120, 194), Point(119, 199), Point(120, 202), Point(123, 207), Point(127, 211), Point(133, 215), Point(142, 219), Point(148, 220), Point(151, 221) };
		std::vector<Point> rightCB = { Point(117, 132), Point(115, 132), Point(115, 129), Point(117, 129), Point(119, 128), Point(122, 127), Point(125, 127), Point(127, 127), Point(130, 127), Point(133, 129), Point(136, 129), Point(138, 130), Point(140, 131), Point(143, 134), Point(144, 136), Point(145, 139), Point(145, 142), Point(145, 145), Point(145, 147), Point(145, 149), Point(144, 152), Point(142, 157), Point(141, 160), Point(139, 163), Point(137, 166), Point(135, 167), Point(133, 169), Point(131, 172), Point(128, 173), Point(126, 176), Point(125, 178), Point(125, 180), Point(125, 182), Point(126, 184), Point(128, 187), Point(130, 187), Point(132, 188), Point(135, 189), Point(140, 189), Point(145, 189), Point(150, 187), Point(155, 186), Point(157, 185), Point(159, 184), Point(156, 185), Point(154, 185), Point(149, 185), Point(145, 187), Point(141, 188), Point(136, 191), Point(134, 191), Point(131, 192), Point(129, 193), Point(129, 195), Point(129, 197), Point(131, 200), Point(133, 202), Point(136, 206), Point(139, 211), Point(142, 215), Point(145, 220), Point(147, 225), Point(148, 231), Point(147, 239), Point(144, 244), Point(139, 248), Point(134, 250), Point(126, 253), Point(119, 253), Point(115, 253) };
		std::vector<Point> star = { Point(75, 250), Point(75, 247), Point(77, 244), Point(78, 242), Point(79, 239), Point(80, 237), Point(82, 234), Point(82, 232), Point(84, 229), Point(85, 225), Point(87, 222), Point(88, 219), Point(89, 216), Point(91, 212), Point(92, 208), Point(94, 204), Point(95, 201), Point(96, 196), Point(97, 194), Point(98, 191), Point(100, 185), Point(102, 178), Point(104, 173), Point(104, 171), Point(105, 164), Point(106, 158), Point(107, 156), Point(107, 152), Point(108, 145), Point(109, 141), Point(110, 139), Point(112, 133), Point(113, 131), Point(116, 127), Point(117, 125), Point(119, 122), Point(121, 121), Point(123, 120), Point(125, 122), Point(125, 125), Point(127, 130), Point(128, 133), Point(131, 143), Point(136, 153), Point(140, 163), Point(144, 172), Point(145, 175), Point(151, 189), Point(156, 201), Point(161, 213), Point(166, 225), Point(169, 233), Point(171, 236), Point(174, 243), Point(177, 247), Point(178, 249), Point(179, 251), Point(180, 253), Point(180, 255), Point(179, 257), Point(177, 257), Point(174, 255), Point(169, 250), Point(164, 247), Point(160, 245), Point(149, 238), Point(138, 230), Point(127, 221), Point(124, 220), Point(112, 212), Point(110, 210), Point(96, 201), Point(84, 195), Point(74, 190), Point(64, 182), Point(55, 175), Point(51, 172), Point(49, 170), Point(51, 169), Point(56, 169), Point(66, 169), Point(78, 168), Point(92, 166), Point(107, 164), Point(123, 161), Point(140, 162), Point(156, 162), Point(171, 160), Point(173, 160), Point(186, 160), Point(195, 160), Point(198, 161), Point(203, 163), Point(208, 163), Point(206, 164), Point(200, 167), Point(187, 172), Point(174, 179), Point(172, 181), Point(153, 192), Point(137, 201), Point(123, 211), Point(112, 220), Point(99, 229), Point(90, 237), Point(80, 244), Point(73, 250), Point(69, 254), Point(69, 252)
		};
		std::vector<Point> pigtail = { Point(81, 219), Point(84, 218), Point(86, 220), Point(88, 220), Point(90, 220), Point(92, 219), Point(95, 220), Point(97, 219), Point(99, 220), Point(102, 218), Point(105, 217), Point(107, 216), Point(110, 216), Point(113, 214), Point(116, 212), Point(118, 210), Point(121, 208), Point(124, 205), Point(126, 202), Point(129, 199), Point(132, 196), Point(136, 191), Point(139, 187), Point(142, 182), Point(144, 179), Point(146, 174), Point(148, 170), Point(149, 168), Point(151, 162), Point(152, 160), Point(152, 157), Point(152, 155), Point(152, 151), Point(152, 149), Point(152, 146), Point(149, 142), Point(148, 139), Point(145, 137), Point(141, 135), Point(139, 135), Point(134, 136), Point(130, 140), Point(128, 142), Point(126, 145), Point(122, 150), Point(119, 158), Point(117, 163), Point(115, 170), Point(114, 175), Point(117, 184), Point(120, 190), Point(125, 199), Point(129, 203), Point(133, 208), Point(138, 213), Point(145, 215), Point(155, 218), Point(164, 219), Point(166, 219), Point(177, 219), Point(182, 218), Point(192, 216), Point(196, 213), Point(199, 212), Point(201, 211)
		};

		std::vector<std::string> names = { "triangle","x","rectangle","circle","check","caret","zig zag","arrow","left square bracket",
			"right square bracket","v","delete","left curly brace","right curly brace","star","pigtail" };
		//Pushing the points with Unistroke initialization 
		this->Unistrokes.push_back(Unistroke(names[0], trianglePoints));
		this->Unistrokes.push_back(Unistroke(names[1], xPoints));
		this->Unistrokes.push_back(Unistroke(names[2], xrectanglePoints));
		this->Unistrokes.push_back(Unistroke(names[3], circle));
		this->Unistrokes.push_back(Unistroke(names[4], checkPoints));
		this->Unistrokes.push_back(Unistroke(names[5], caret));
		this->Unistrokes.push_back(Unistroke(names[6], zigZag));
		this->Unistrokes.push_back(Unistroke(names[7], arrow));
		this->Unistrokes.push_back(Unistroke(names[8], leftsquarebracket));
		this->Unistrokes.push_back(Unistroke(names[9], rightsquarebracket));
		this->Unistrokes.push_back(Unistroke(names[10], v));
		this->Unistrokes.push_back(Unistroke(names[11], deletePoints));
		this->Unistrokes.push_back(Unistroke(names[12], leftCB));
		this->Unistrokes.push_back(Unistroke(names[13], rightCB));
		this->Unistrokes.push_back(Unistroke(names[14], star));
		this->Unistrokes.push_back(Unistroke(names[15], pigtail));

	}
	std::vector<Unistroke> Unistrokes;
	std::vector<Unistroke> OfflineStrokes;

	GestureRecognizer(vector<Unistroke>& inputTempatePoints) {

		this->OfflineStrokes = inputTempatePoints;

	}
	Result Recognize(vector<Point>& points, bool useProtractor) {
		//TODO Raghav Vikranth
		auto t0 = std::chrono::high_resolution_clock::now();
		string str = "";
		Unistroke candidate(str, points);//resampled here for candiate or input gesture
		vector<double> nbest;
		int u = -1;//If nothing matches the initializaiton
		double b = INFINITY; //intMAX
		for (int i = 0; i < Unistrokes.size(); i++) {
			double d;

			if (useProtractor)//this is never invoked for now.
				d = OptimalCosineDistance(Unistrokes[i].vectorizedPoints, candidate.vectorizedPoints);
			else
				d = DistanceAtBestAngle(candidate.points, Unistrokes[i], -AngleRange, AngleRange, AnglePrecision);

			double tempscore = (1.0 - d / HalfDiagonal);
			nbest.push_back(tempscore);
			if (d < b) {
				b = d;
				u = i; //its going to point at the template that is under consideration and saves the value for the same for the min distance
			}
		}
		sort(nbest.begin(), nbest.end());
		auto t1 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();//Time elapsed
		return (u == -1) ? Result("No match.", 0.0, duration) :
			Result(Unistrokes[u].name, useProtractor ? (1.0 - b) : (1.0 - b / HalfDiagonal), duration);
	}
	bool comparePair(const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) {
		return a.second > b.second;
	}

	OfflineResult OfflineRecognize(vector<Point>& candidatePoints, bool useProtactor) {
		auto t0 = std::chrono::high_resolution_clock::now();
		string str = "";
		//Unistroke candidate(str, candidatePoints);//resampled here for candiate or input gesture
		vector<OfflineResult> final;
		int u = -1;//If nothing matches the initializaiton
		double b = INFINITY; //intMAX
		OfflineResult temp;
		vector<pair<string, double>> t1;

		for (int i = 0; i < OfflineStrokes.size(); i++) {
			double d;
			d = DistanceAtBestAngle(candidatePoints, OfflineStrokes[i], -AngleRange, AngleRange, AnglePrecision);

			temp.gestureName = OfflineStrokes[i].name;
			temp.score = 1.0 - d / HalfDiagonal;
			string name = temp.gestureName;
			t1.push_back(make_pair(name, temp.score));
			//need to sort TODO
			if (d < b) {
				b = d;
				u = i; //its going to point at the template that is under consideration and saves the value for the same for the min distance
			}
		}
		temp.gestureName = OfflineStrokes[u].name;
		temp.score = (1.0 - b / HalfDiagonal);
		std::sort(t1.begin(), t1.end(), std::bind(&GestureRecognizer::comparePair, this, std::placeholders::_1, std::placeholders::_2));
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
			vector<string> labelList = { "triangle","x","rectangle","circle","check","caret","arrow","left_sq_bracket",
				"right_sq_bracket","v","delete_mark","left_curly_brace","right_curly_brace","star","pigtail","question_mark" };
			string fileName = "";
			string part = "";
			string storeName = "";
			for (int k = 2; k <= 11; k++) {
				if (k < 10) {
					part = "C:/Users/vikra/Desktop/HCIRA/xml/xml_logs/s0" + to_string(k);
					storeName = "s0" + to_string(k);
				}
				else {
					part = "C:/Users/vikra/Desktop/HCIRA/xml/xml_logs/s" + to_string(k);
					storeName = "s" + to_string(k);
				}

				for (int j = 0; j < 16; j++) { //this is for the gestures
					vector<vector<Point>> listOfPoints;
					const char* gestureName; //initialized the gesture name
					for (int i = 1; i <= 10; i++) {
						vector<Point> points;
						if (i < 10) {
							fileName = part + "/medium/" + labelList[j] + "0" + to_string(i) + ".xml";
						}
						else {
							fileName = part + "/medium/" + labelList[j] + to_string(i) + ".xml";
						}

						doc.LoadFile(fileName.c_str()); //C:\Users\91977\Downloads\xml\xml_logs\s02\medium\arrow01.xml
						//wxLogMessage(fileName.c_str());
						// Get the root element of the XML document
						tinyxml2::XMLElement* root = doc.FirstChildElement("Gesture");


						// Get the value of attributes in the root element
						gestureName = root->Attribute("Name");
						const char* appName = root->Attribute("AppName");
						int numPts = root->IntAttribute("NumPts");
						//wxLogMessage(gestureName);
						//wxLogMessage(appName);

						//wxLogMessage(numPts);

						//std::cout << "Gesture name: " << gestureName << std::endl;
						//std::cout << "Application name: " << appName << std::endl;
						//std::cout << "Number of points: " << numPts << std::endl;
						 //pushing points into this from xml

						// Iterate through the child elements (i.e., the points) of the root element
						for (tinyxml2::XMLElement* point = root->FirstChildElement(); point != NULL; point = point->NextSiblingElement())
						{
							// Get the values of the X, Y, and T attributes of the point element
							double x = point->IntAttribute("X");
							double y = point->IntAttribute("Y");
							double t = point->IntAttribute("T");
							/*wxLogMessage("x co ordinate %f", x);
							wxLogMessage("y co ordinate %f", y);
							wxLogMessage("Time %f", t);*/
							Point tempPoint(x, y);
							points.push_back(tempPoint);
							//std::cout << "X: " << x << ", Y: " << y << ", T: " << t << std::endl;
						}
						listOfPoints.push_back(points);


						//offlineData[""][""][""].push_back({ tempPoint });
						//1.Capture the points from XML 
						//2. Resample them with with unistroke and add gesture label 
						//3. Store these points in the the class instance of GestureRecognizer
						//offlineData[storeName]["medium"][labelList[j]] = listOfPoints;
					}
					offlineData[storeName]["medium"][labelList[j]] = listOfPoints;
				}
				//wxLogMessage("offlineData:");
			}
			//for (const auto& it1 : offlineData) {
			//	wxLogMessage("Key 1: %s", it1.first);

			//	for (const auto& it2 : it1.second) {
			//		wxLogMessage("Key 2: %s", it2.first);

			//		for (const auto& it3 : it2.second) {
			//			wxLogMessage("Key 3: %s", it3.first);

			//			for (const auto& vec : it3.second) {
			//				for (const auto& v : vec) {
			//					// Loop through each Point object in the current vector
			//						// Print the x and y coordinates of the current Point object to the log
			//						wxLogMessage("Point: (%f, %f)", v.x, v.y);
			//					
			//				}
			//			}
			//		}
			//	}
			//}
		}


	void preProcessOfflineData() {
		preProcessedData = offlineData;

		for (auto& user : offlineData) {
			for (auto& speed : user.second) {
				if (speed.first == "medium") {
					for (auto& gesture : speed.second) {
						//wxLogMessage("");
						preProcessedData[user.first][speed.first][gesture.first] = vector< vector<Point>>();
						for (auto& temp : gesture.second) {
							//gesture name , points
							string s = gesture.first;
							Unistroke unistroke(s, temp);
							preProcessedData[user.first][speed.first][gesture.first].push_back({ unistroke.points });
						}
					}
				}
			}
		}
	}


	void recognizeOfflineData() {

		map< string, map<int, map< string, double>>> score;
		double total = 0;
		double correct = 0;
		std::ofstream outfile("output.csv");
		LogData log;

		if (!outfile.is_open()) {
			std::cerr << "Failed to open file for writing." << std::endl;
		}

		outfile << "User[all-users],GestureType[all-gesture-types],RandomIteration[1to10],#ofTrainingExample[E],TotalSizeOfTrainingSet[count],TrainingSetContents[specific-gesture-instances],Candidate[specific-instance],RecoResultGestureType[what-was-recognized],CorrectIncorrect[1or0],RecoResultScore,RecoResultBestMatch[specific-instance],RecoResultNBestSorted[instance-and-score]" << std::endl;
		for (auto& user : preProcessedData) {
		
			{
				log.User = user.first;
				for (int gesture = 1; gesture <= 9; gesture++) { //this is the value of the 

					for (int i = 1; i <= 10; i++) { //we can decrease the value to 10

						log.RandomIteration = to_string(i);
						log.NoTrainingExample = gesture;
						pair<vector<Unistroke>, map< string, vector<Point>>> split_data = getSplitData(preProcessedData[user.first]["medium"], gesture);

						vector<Unistroke> train_set = split_data.first;//TEsting conatines everything from offline 160-training
						map<string, vector<Point>> test_set = split_data.second;//a specific set of gestures with points training


						log.TotalSizeOfTrainingSet = train_set.size();
						log.TrainingSetContents.clear();
						for (auto& elem : train_set) {
							string temp = user.first + "-" + elem.name;
							log.TrainingSetContents.push_back(temp);
						}


						GestureRecognizer recognizer(train_set);//loading templates

						
						for (const auto& elem : test_set) {
							//wxLogMessage("Gesture name in test set %s",elem.first);
							
							std::string r1 = ""; //gesturType
							
							for (char c : elem.first) {
								if (!isdigit(c)) {
									r1 += c;
								}
							}
							
							
							vector<Point> pts = elem.second;//candidates pts
							OfflineResult res = recognizer.OfflineRecognize(pts, false);

							string gestureRecognized = res.gestureName;
							vector<pair<string, double>> Nbest = res.nbest;

							//
							std::string r2 = ""; //Recoresult gesture type

							for (char c : gestureRecognized) {
								if (!isdigit(c)) {
									r2 += c;
								}
							}
							/*r1.erase(r1.begin(), std::find_if(r1.begin(), r1.end(), [](unsigned char ch) {
								return !std::isspace(ch);
								}));
							r1.erase(std::find_if(r1.rbegin(), r1.rend(), [](unsigned char ch) {
								return !std::isspace(ch);
								}).base(), r1.end());
							r2.erase(r2.begin(), std::find_if(r2.begin(), r2.end(), [](unsigned char ch) {
								return !std::isspace(ch);
								}));
							r2.erase(std::find_if(r2.rbegin(), r2.rend(), [](unsigned char ch) {
								return !std::isspace(ch);
								}).base(), r2.end());*/

							log.GestureType = r1;
							log.RecoResultGestureType = r2;
							//wxLogMessage("%s", r2);
							if (log.GestureType == log.RecoResultGestureType) {
								log.correct = true ;
							}
							else {
								log.correct = false;
							}
							//log.correct = gestureRecognized.substr(0, elem.first.length() - 1) == elem.first.substr(0, elem.first.length() - 1);
							log.RecoResultScore = to_string(res.score);
							log.RecoResultBestMatch = user.first + "-" + res.gestureName;
							log.RecoResultNBest = res.nbest;
							log.candidateSpecificInstance = user.first + "-" + elem.first;

							if (log.GestureType==log.RecoResultGestureType) {

								score[user.first][gesture][elem.first] += 1;

								//log.GestureType = elem.first;
								correct += 1.0;

							}

							total += 1.0;

							
							std::string contentsString = vector1ToString(log.TrainingSetContents);

							std::string nbestStr = vectorToString(log.RecoResultNBest);
							outfile << log.User << "," << log.GestureType << "," << log.RandomIteration << "," << log.NoTrainingExample << "," << log.TotalSizeOfTrainingSet << "," << contentsString << "," << log.candidateSpecificInstance << "," << log.RecoResultGestureType << "," << to_string(log.correct) << "," << log.RecoResultScore << "," << log.RecoResultBestMatch << "," << nbestStr << endl;

						}

					}

				}
			}

		}
		wxLogMessage("total %f", total);
		wxLogMessage("correct %f", correct);

		outfile << "The Total Accuracy" << endl;
		outfile << to_string((correct / total) * 100.0) << endl;

		outfile.close();


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

	
	//map< string, vector< vector<Point>>>
	pair<vector<Unistroke>, map< string, vector<Point>>> getSplitData(map< string, vector< vector<Point>>>& gestures, int E) {
		//map< string, vector< vector<Point>>> training_set;
		vector<Unistroke> training_set;
		map< string, vector<Point>> testing_set;
		srand(time(nullptr));
		//we select the testing set here
		//randomly 
		//and then place if an condition in the for loop when we hit the random value we continue and do nothing
		for (auto& gesture : gestures) {
			for (int i = 1; i <= E; i++) {
				string name = gesture.first + to_string(i);
				training_set.push_back(Unistroke(name, gesture.second[i - 1]));
				//training_set[gesture.first] = temp;
			}
			int range = 10 - E;
			int random_val=rand() % range + E ;
			//int s = rand() % (10 - E) + E+1;
			string name = gesture.first + to_string(random_val);
			testing_set[name] = gesture.second[random_val];
		}
		return  make_pair(training_set, testing_set);
	}

};






class MyCanvas : public wxWindow
{
public:
	//Canvas for drawing 
	MyCanvas(wxWindow* parent) : wxWindow(parent, wxID_ANY)
	{
		Connect(wxEVT_PAINT, wxPaintEventHandler(MyCanvas::OnPaint));
		Connect(wxEVT_LEFT_DOWN, wxMouseEventHandler(MyCanvas::OnLeftDown));
		Connect(wxEVT_LEFT_UP, wxMouseEventHandler(MyCanvas::OnLeftUp));
		Connect(wxEVT_MOTION, wxMouseEventHandler(MyCanvas::OnMotion));
		m_output = new wxStaticText(this, wxID_ANY, "");
		wxFont font(20, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
		m_font = font;
		m_output->SetFont(m_font);
		//XML Files points to be loaded
		OfflineRecognizer a;
		a.preProcessOfflineData();
		a.recognizeOfflineData();

	}
private:
	//event handler is called when the canvas is repainted
	void OnPaint(wxPaintEvent& event)
	{
		wxPaintDC dc(this);
		wxGraphicsContext* gc = wxGraphicsContext::Create(dc);

		gc->SetPen(wxPen(*wxRED, 5));
		for (size_t i = 1; i < m_points.size(); i++)
		{
			dc.DrawLine(m_points[i - 1].x, m_points[i - 1].y, m_points[i].x, m_points[i].y);
		}

		delete gc;
	}
	//event handler is called when the left mouse button is pressed and pushes position to m_points
	void OnLeftDown(wxMouseEvent& event)
	{
		m_points.clear();
		m_points.push_back(event.GetPosition());

		//a.Recognize(m_points,true); did't use this again 
		CaptureMouse();
	}
	//event handler is called when the left mouse button is released
	void OnLeftUp(wxMouseEvent& event)
	{
		m_points.push_back(event.GetPosition());
		ReleaseMouse();
		points.clear();
		for (auto& wxPoint : m_points) {
			points.emplace_back(wxPoint.x, wxPoint.y);
		}
		GestureRecognizer GR;//Our $1 recognizer object


		Result res = GR.Recognize(points, false);
		wxString outputStr = wxString::Format("Result: %s (%f) in %d ms", res.Name, res.Score, res.Time);
		m_output->SetLabel(outputStr);


		//for (auto point : points) {
		//    cout << "x: " << point.x << " y: " << point.y << endl; //this didn't write the output on the console rather has to set a breakpoint to check the values of the variables.
		//}
		Refresh();
	}
	//function event handler is called when the mouse is moved
	void OnMotion(wxMouseEvent& event)
	{
		if (event.Dragging() && event.LeftIsDown())
		{
			m_output->SetLabel("Recording unistroke...");
			m_points.push_back(event.GetPosition());
			Refresh();
		}
	}

	vector<wxPoint> m_points;
	//points being stored for future reference
	vector<Point> points;
	wxStaticText* m_output;

};



//
//class MyApp : public wxApp
//{
//public:
//	virtual bool OnInit()
//	{
//		wxFrame* frame = new wxFrame(NULL, wxID_ANY, "$1 Recognizer");
//		MyCanvas* canvas = new MyCanvas(frame);
//		frame->Show();
//
//
//		return true;
//	}
//
//	void LoadVal() {
//		//Importing points from XML to my data set of points
//		//16 gestures from 10 users 
//		//Total 11 files
//
//
//	}
//
//
//
//
//
//
//
//};
//
//
////creates an instance of the MyApp
//wxIMPLEMENT_APP(MyApp);
