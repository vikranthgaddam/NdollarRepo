#include <wx/wx.h>
#include <wx/graphics.h>]
#include <vector>
#include <random>

using namespace std;

//class definitions
//Defining point struct
//class Point { 
//public:
//	double x;
//	double y;
//	Point(double x, double y) {
//		this->x = x;
//		this->y = y;
//	}
//};

//class MyCanvas : public wxWindow
//{
//public:
//	//Canvas for drawing 
//	MyCanvas(wxWindow* parent) : wxWindow(parent, wxID_ANY)
//	{
//
//		m_prompt = new wxStaticText(this, wxID_ANY, "Draw a {Gesture}");
//		//m_output = new wxStaticText(this, wxID_ANY, "");
//		total_counter = new wxStaticText(this, wxID_ANY, "Total Counter: 0/160");
//		
//		m_clearButton = new wxButton(this, wxID_ANY, "Clear");
//		m_submitButton = new wxButton(this, wxID_ANY, "Submit");
//
//		wxFont font(10, wxFONTFAMILY_DEFAULT, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL);
//		m_font = font;
//		//m_output->SetFont(m_font);
//		
//		total_counter->SetFont(m_font);
//		m_prompt->SetFont(m_font);
//
//		m_submitButton->SetBackgroundColour(wxColour(0, 255, 0));
//		m_clearButton->SetBackgroundColour(wxColour(255, 0, 0));
//
//		//Main Sizer 
//		wxBoxSizer* sizer = new wxBoxSizer(wxVERTICAL);
//		
//		//sizer->Add(m_output, 0, wxALIGN_CENTER_HORIZONTAL | wxTOP, 0);
//		//Text Sizer
//		wxBoxSizer* textSizer = new wxBoxSizer(wxHORIZONTAL);
//		textSizer->Add(m_prompt, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 0);
//		textSizer->Add(total_counter, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 0); // add counter to sizer
//
//		//Button Sizer
//		wxBoxSizer* buttonsSizer = new wxBoxSizer(wxHORIZONTAL);
//		buttonsSizer->Add(m_clearButton, 0, wxALIGN_CENTER_VERTICAL | wxRIGHT, 20);
//		buttonsSizer->Add(m_submitButton, 0, wxALIGN_CENTER_VERTICAL | wxLEFT, 20);
//		//creating space from text
//		sizer->Add(textSizer, 0, wxALIGN_CENTRE_HORIZONTAL | wxTop, 10);
//		sizer->AddStretchSpacer(1);
//		sizer->Add(buttonsSizer, 0, wxALIGN_CENTRE_HORIZONTAL | wxBottom, 10);
//
//
//		
//		//sizer->Add(m_clearButton, 0, wxALIGN_CENTER_HORIZONTAL | wxBottom, 0);
//		//sizer->Add(m_submitButton, 0, wxALIGN_CENTER_HORIZONTAL | wxBottom, 0);
//		SetSizerAndFit(sizer);
//
//		Connect(wxEVT_PAINT, wxPaintEventHandler(MyCanvas::OnPaint));
//		Connect(wxEVT_LEFT_DOWN, wxMouseEventHandler(MyCanvas::OnLeftDown));
//		Connect(wxEVT_LEFT_UP, wxMouseEventHandler(MyCanvas::OnLeftUp));
//		Connect(wxEVT_MOTION, wxMouseEventHandler(MyCanvas::OnMotion));
//		Connect(m_clearButton->GetId(), wxEVT_BUTTON, wxCommandEventHandler(MyCanvas::OnClearButtonClick));
//		Connect(m_submitButton->GetId(), wxEVT_BUTTON, wxCommandEventHandler(MyCanvas::OnSubmitButtonClick));
//
//	}
//private:
//	wxStaticText* m_prompt;
//	wxStaticText* total_counter;
//	wxButton* m_clearButton;
//	wxButton* m_submitButton;
//	vector<string> labelList = { "triangle","x","rectangle","circle","check","caret","arrow","left_sq_bracket",
//			"right_sq_bracket","v","delete_mark","left_curly_brace","right_curly_brace","star","pigtail","zig-zag" };
//	unordered_map<string, int> mp;
//	int counter = 0;
//	int current_gesture = 0;
//	
//	//event handler is called when the canvas is repainted
//	void OnPaint(wxPaintEvent& event)
//	{
//		wxPaintDC dc(this);
//		wxGraphicsContext* gc = wxGraphicsContext::Create(dc);
//
//		gc->SetPen(wxPen(*wxRED, 5));
//		for (size_t i = 1; i < m_points.size(); i++)
//		{
//			dc.DrawLine(m_points[i - 1].x, m_points[i - 1].y, m_points[i].x, m_points[i].y);
//		}
//
//		delete gc;
//	}
//	// event handler is called when the "Clear" button is clicked
//	void OnClearButtonClick(wxCommandEvent& event)
//	{
//		m_points.clear();
//		Refresh();
//	}
//
//	void OnSubmitButtonClick(wxCommandEvent& event) {
//		m_points.clear();
//		Refresh();
//		wxString counterStr = wxString::Format("Counter: %d", counter);
//		total_counter->SetLabel(counterStr);
//	}
//
//	//event handler is called when the left mouse button is pressed and pushes position to m_points
//	void OnLeftDown(wxMouseEvent& event)
//	{
//		m_points.push_back(event.GetPosition());
//		CaptureMouse();
//	}
//	//event handler is called when the left mouse button is released
//	void OnLeftUp(wxMouseEvent& event)
//	{
//		m_points.push_back(event.GetPosition());
//		points.clear();
//		ReleaseMouse();
//		Refresh();
//	}
//	//function event handler is called when the mouse is moved
//	void OnMotion(wxMouseEvent& event)
//	{
//		if (event.Dragging() && event.LeftIsDown())
//		{
//			m_points.push_back(event.GetPosition());
//			Refresh();
//		}
//	}
//
//	vector<wxPoint> m_points;
//	//points being stored for future reference
//	vector<Point> points;
//	wxStaticText* m_output;
//
//};
//
//
//
//
