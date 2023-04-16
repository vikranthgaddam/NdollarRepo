//Author Raghav Gupta and Vikranth Gaddam
//CPP of $1 Recognizer
#include <wx/wx.h>
#include <wx/graphics.h>
#include <iostream>
#include <tinyxml2.h>


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

class MyCanvas : public wxWindow
{
public:
	//Canvas for drawing 
	MyCanvas(wxWindow* parent) : wxWindow(parent, wxID_ANY)
	{

		m_prompt = new wxStaticText(this, wxID_ANY, "Draw a {Gesture}");
		//m_output = new wxStaticText(this, wxID_ANY, "");
		total_counter = new wxStaticText(this, wxID_ANY, "Total Counter: 0/160");
		SetBackgroundStyle(wxBG_STYLE_PAINT);
		m_clearButton = new wxButton(this, wxID_ANY, "Clear");
		m_submitButton = new wxButton(this, wxID_ANY, "Submit");

		wxFont font(10, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
		m_font = font;
		//m_output->SetFont(m_font);

		total_counter->SetFont(m_font);
		m_prompt->SetFont(m_font);

		m_submitButton->SetBackgroundColour(wxColour(0, 255, 0));
		m_clearButton->SetBackgroundColour(wxColour(255, 0, 0));

		//Main Sizer 
		wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);

		//sizer->Add(m_output, 0, wxALIGN_CENTER_HORIZONTAL | wxTOP, 0);
		//Text Sizer
		wxBoxSizer* textSizer = new wxBoxSizer(wxHORIZONTAL);
		textSizer->Add(m_prompt, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 0);
		textSizer->Add(total_counter, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 0); // add counter to sizer

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
		
	}
private:
	
	wxStaticText* m_prompt;
	wxStaticText* total_counter;
	wxButton* m_clearButton;
	bool m_isDrawing = false;
	wxButton* m_submitButton;
	vector<wxPoint> m_points;
	vector<vector<wxPoint>> listOf;
	vector<vector<Point>> strokes;
	wxFont m_font;
	vector<Point> points;
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
			/*wxLogMessage("Strokes:");
			for (size_t i = 0; i < strokes.size(); i++) {
				wxLogMessage("Stroke %d:", i);
				const vector<Point>& points = strokes[i];
				for (size_t j = 0; j < points.size(); j++) {
					wxLogMessage("    (%f, %f)", points[j].x, points[j].y);
				}
			}*/
			
			m_points.clear();
			Refresh();
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
			//m_output->SetLabel("Recording unistroke...");
			m_points.push_back(event.GetPosition());
			Refresh();
		}
	}
	

};

class MyApp : public wxApp
{
public:
	virtual bool OnInit();
};

bool MyApp::OnInit(){

	wxFrame* frame = new wxFrame(NULL, wxID_ANY, "$1 Recognizer");
	MyCanvas* canvas = new MyCanvas(frame);
	return true;
}

wxIMPLEMENT_APP(MyApp);