#include "cash-display.hpp"

extern int scale;

bool CashDisplay::no_other_instantiation = true;


CashDisplay::CashDisplay(const int _window_row,const int _window_col, const std::vector<CashPanelInfo>& _panel_info,const int _scale,const unsigned char _margin_color)
  : window_row(_window_row),
    window_col(_window_col),
    window_is_open(false),
    png_is_open(false),
    margin_color(_margin_color)
{
  /* CashDisplay must have only one instantiation at a time */
  if(no_other_instantiation){
    no_other_instantiation = false;
  }else{
    std::cerr << "CashDisplay::CashDisplay() the number of instantiated Display objects cannot be more than two at the same time." << std::endl;
    exit(-1);
  }

  /* Allocate memory in data and initialize it. */
  if(window_row==0 || window_col==0){
    std::cerr << "CashDisplay::Data::Data() Error: nrow=ncol=0 is not allowed." << std::endl;
    return;
  }
  try{
    data = new unsigned char [window_row*window_col];
  }catch(std::bad_alloc){
    std::cerr << "CashDisplay::Data::Data() Memory exausted" << std::endl;
    exit(-1);
  }
  memset(data,margin_color,sizeof(unsigned char)*window_row*window_col);

  /* An error check */
  if(_panel_info.empty()){
    std::cerr << "CashDisplay::CashDisplay() the number of elements in panel_info is zero." << std::endl;
    exit(-1);
  }

  /* Get each panel information */
  for(std::vector<CashPanelInfo>::const_iterator iter=_panel_info.begin();iter!=_panel_info.end();iter++){
    if(iter->n_row<=0 ||
       iter->n_col<=0 ||
       iter->o_row<0  ||
       iter->o_col<0){
      std::cerr << "CashDisplay::CashDisplay() the " << _panel_info.size()-static_cast<std::vector<CashPanelInfo>::size_type>(_panel_info.end()-iter) << "-th panel_info element has abnormal specification: n_row=" << iter->n_row << " n_col=" << iter->n_col << " o_row=" << iter->o_row << " o_col=" << iter->o_col << std::endl;
      exit(-1);
    }
   
    if(iter->n_row+iter->o_row > window_row ||
       iter->n_col+iter->o_col > window_col){
      std::cerr << "CashDisplay::CashDisplay() the " << _panel_info.size()-static_cast<std::vector<CashPanelInfo>::size_type>(_panel_info.end()-iter) << "-th panel_info element's specification violates the window size: n_row=" << iter->n_row << " n_col=" << iter->n_col << " o_row=" << iter->o_row << " o_col=" << iter->o_col << " window_row=" << window_row << " window_col=" << window_col << std::endl;
      exit(-1);
    }

    panels.push_back(*iter);
  }

  /* Set the size of pixel */
  if(_scale>0){
    scale = _scale;
  }
  else{
    std::cerr << "CashDisplay::CashDisplay(): Error, scale=" << _scale << std::endl;
    exit(-1);
  }

  /* Set default colors */
  ColorRGB(0, 0x0, 0x0, 0x0); // dead
  ColorRGB(1, 0x1a,0x95,0xcc); // bacteria - BLUE
  ColorRGB(2, 0x11,0x58,0x79); // bacteria, plasmid - BLUE DARK
  ColorRGB(3, 0x47,0xba,0xfb); // bacteria, resistant plasmid - BLUE LIGHT

  ColorRGB(4, 0xdd,0x5a,0x3a); // resistant bacteria - RED
  ColorRGB(5, 0x88,0x3a,0x1a); // resistant bacteria, plasmid - RED DARK
  ColorRGB(6, 0xff,0x7a,0x5a); // resistant bacteria, resistant plasmid - RED LIGHT

  ColorRGB(7, 0x6b,0x2f,0xd9); // bacteria, phage - GREENBLUE
  ColorRGB(8, 0x41,0x1a,0x83); // bacteria, plasmid, phage - GREENBLUE DARK
  ColorRGB(9, 0x9a,0x47,0xff);  // bacteria, resistant plasmid, phage - GREENBLUE LIGHT 

  ColorRGB(10, 0x1a,0xbe,0x6b); // resistant bacteria, phage - PURPLE
  ColorRGB(11, 0x11,0x71,0x3e); // resistant bacteria, plasmid, phage - PURPLE DARK
  ColorRGB(12, 0x3e,0xf5,0x9d); // resistant bacteria, resistant plasmid, phage - PURPLE LIGHT 

  ColorRGB(13, 0xa4,0x39,0x8c); // bacteria, resistant phage - MAROON
  ColorRGB(14, 0x62,0x21,0x55); // bacteria, plasmid, resistant phage - MAROON DARK
  ColorRGB(15, 0xdc,0x6e,0xc0); // bacteria, resistant plasmid, resistant phage - MAROON LIGHT


  ColorRGB(16, 0xc8,0xce,0x00); // resistant bacteria, resistant phage - YELLOW
  ColorRGB(17, 0x76,0x7f,0x00); // resistant bacteria, plasmid, resistant phage - YELLOW GREEN DARK
  ColorRGB(18, 0xff,0xff,0x4c); // resistant bacteria, resistant plasmid, resistant phage - YELLOW LIGHT

  ColorRGB(19, 0xff,0xff,0xff); // virion
  ColorRGB(20, 0xff,0xff,0xff); // resistant virion
  ColorRGB(21, 0x88,0x88,0x88); // exposed cell 0x105
  ColorRGB(22, 0xFF,0x88,0x88); // exposed cell 0x115
  ColorRGB(23, 0x88,0xFF,0x88); // exposed cell 0x195
  ColorRGB(24, 0x88,0x88,0xFF); // exposed cell 0x905
  ColorRGB(25, 0xFF,0xAA,0x88); // exposed cell 0x915
  ColorRGB(26, 0xFF,0xFF,0x88); // exposed cell 0x995
  ColorRGB(27, 0xFF,0x88,0xAA); // exposed cell 0x10D
  ColorRGB(28, 0xFF,0x88,0xFF); // exposed cell 0x11D
  ColorRGB(29, 0x88,0xFF,0xAA); // exposed cell 0x19D
  ColorRGB(30, 0x88,0xFF,0xFF); // exposed cell 0x90D
  ColorRGB(31, 0xAA,0xFF,0x88); // exposed cell 0x91D
  ColorRGB(32, 0x00,0xFF,0x88); // exposed cell 0x99D
  // ColorRGB(43, 0xFF,0xD5,0X00); //high virion concentration
  ColorRGB(43, 0xFF,0xD5,0X00); //high virion concentration
  // ColorRGB(45, 0xDD,0xB3,0x00); //moderate high virion concentration
  // ColorRGB(42, 0xBB,0x8C,0x00); //moderate high virion concentration
  ColorRGB(42, 0xaa,0x76,0x00); //moderate high virion oconcentration
  // ColorRGB(43, 0x99,0x66,0x00); //moderate virion concentration
  ColorRGB(41, 0x66,0x3c,0x00); //moderate low virion concentration
  // ColorRGB(41, 0x55,0x34,0x00); //low virion concentration
  ColorRGB(40, 0x00,0x00,0x00); //no virions
  
  // ColorRGB(CashColor::BLACK,0, 0, 0);
  // ColorRGB(CashColor::WHITE,255, 255, 255);
  // ColorRGB(CashColor::RED,255, 0, 0);
  // ColorRGB(CashColor::GREEN,0, 255, 0);
  // ColorRGB(CashColor::BLUE,0, 0, 255);
  // ColorRGB(CashColor::YELLOW,255, 255, 0);
  // ColorRGB(CashColor::BROWN,188, 143, 143);
  // ColorRGB(CashColor::GRAY,220, 220, 220);
  // ColorRGB(CashColor::VIOLET,148, 0, 211);
  // ColorRGB(CashColor::CYAN,0, 255, 255);
  // ColorRGB(CashColor::MAGENTA,255, 0, 255);
  // ColorRGB(CashColor::ORANGE,255, 165, 0);
  // ColorRGB(CashColor::INDIGO,114, 33, 188);
  // ColorRGB(CashColor::MAROON,103, 7, 72);
  // ColorRGB(CashColor::TURQUOISE,64, 224, 208);
  // ColorRGB(CashColor::GREEN4,0, 139, 0);
}

CashDisplay::~CashDisplay()
{
  delete[] data;

  if(window_is_open)
    CloseDisplayImmediately();
  if(png_is_open)
    ClosePNG();

  no_other_instantiation = true;
}

void CashDisplay::color_yellow2yellow(const unsigned char color_ind_begin,const unsigned char length)
{
  int r=256,g=255,b=0;
  double x = 1530./static_cast<double>(length-1);
  int y=0;
  for(int i=0;i<1531;++i){
    if(i<=255){
      r = r - 1;
    }else if(i<=255*2){
      b = b + 1;
    }else if(i<=255*3){
      g = g - 1;
    }else if(i<=255*4){
      r = r + 1;
    }else if(i<=255*5){
      b = b - 1;
    }else if(i<=255*6){
      g = g + 1;
    }
    if(i == static_cast<int>(y*x+0.5)){
      ColorRGB(y+color_ind_begin,r,g,b);
      y++;
    }
  }
}

void CashDisplay::color_yellow2red(const unsigned char color_ind_begin,const unsigned char length)
{
  int r=256,g=255,b=0;
  double x = 1275./static_cast<double>(length-1);
  int y=0;
  for(int i=0;i<1276;++i){
    if(i<=255){
      r = r - 1;
    }else if(i<=255*2){
      b = b + 1;
    }else if(i<=255*3){
      g = g - 1;
    }else if(i<=255*4){
      r = r + 1;
    }else if(i<=255*5){
      b = b - 1;
    }
    if(i == static_cast<int>(y*x+0.5)){
      ColorRGB(y+color_ind_begin,r,g,b);
      y++;
    }
  }
}

void CashDisplay::put_histogram(const unsigned panel_ind,const std::vector<double>& frequencies,const double max_y,const unsigned char color_of_blank,const unsigned char color_of_background,const std::vector<unsigned char>& color_of_bar)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||panel_ind < panels.size());
  }catch(GeneralError){
    std::cerr << "CashDisplay::put_histogram(): Error, out-of-bound panel_ind=" << panel_ind << std::endl;
    exit(-1);
  }
  
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(max_y>0.));
  }catch(GeneralError){
    std::cerr << "CashDisplay::put_histogram(): Error, max_y=" << max_y << std::endl;
    exit(-1);
  }
  
  int n_pixel_per_bin = panels[panel_ind].n_col/frequencies.size();
  
  if(ASSERT::ERROR_CHECK){
    if(n_pixel_per_bin < 1){
      std::cerr << "CashHistogram::CashHistogram() Error, the number of bins must be smaller than n_row." << std::endl;
      exit(-1);
    }
  }
  
  if(ASSERT::ERROR_CHECK){
    for(std::vector<double>::const_iterator iter=frequencies.begin();iter!=frequencies.end();iter++){
      if(*iter<0.){
	std::cerr << "CashDisplay::put_histogram(): Error, frequencies have negative value." << std::endl;
	exit(-1);
      }
    }
  }

  int row,col,height;
  int col_limit, anti_height;
  std::vector<unsigned char>::const_iterator color_iter = color_of_bar.begin();
  
  /* Draw each bin */
  for(int bin=0; bin<static_cast<int>(frequencies.size()); ++bin){
    height = static_cast<unsigned>(frequencies[bin]/max_y*panels[panel_ind].n_row + 0.5);

    /* If height is out of boundary, set it to maximum */
    if(height > panels[panel_ind].n_row) 
      height = panels[panel_ind].n_row;

    anti_height = panels[panel_ind].n_row - height;

    col_limit = (bin+1)*n_pixel_per_bin+1;
    /* Draw each column */
    for(col=bin*n_pixel_per_bin+1; col<col_limit; ++col){
      /* Draw each row */
      for(row=panels[panel_ind].n_row; row>0; --row){
	data_at(panel_ind,row,col) = (row > anti_height) ? (*color_iter) : color_of_blank;
      }
    }
    /* Incliment the colors. */
    if(color_iter + 1 < color_of_bar.end())
      ++color_iter;
  }

  /* Fill the rest by blank */
  for(col=static_cast<int>(frequencies.size())*n_pixel_per_bin+1; col <= panels[panel_ind].n_col; ++col)
    for(row = panels[panel_ind].n_row; row>0; --row){
      data_at(panel_ind,row,col) = color_of_background;
    }
}

bool CashDisplay::xy_window_to_rc_panel(int x,int y,unsigned& panel_ind,int& row,int& col) const
{
  x = x/scale + 1;
  y = y/scale + 1;

  for(Panels::const_iterator iter=panels.begin();iter!=panels.end();iter++){
    if(iter->o_row < y &&
       iter->o_row +iter->n_row >= y &&
       iter->o_col < x &&
       iter->o_col +iter->n_col >= x){
      panel_ind = static_cast<unsigned>(iter - panels.begin());
      row = y - iter->o_row;
      col = x - iter->o_col;
      return true;
    }
  }
  
  /* if we reach here, there is no panel correspongind to this
     (x,y) coordinate. */
  return false;
}

/******************************************************/

int CashDisplay::get_n_row(const unsigned panel_ind) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||panel_ind < panels.size());
  }catch(GeneralError){
    std::cerr << "CashDisplay::get_n_row(): Error, out-of-bound panel_ind=" << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].n_row;
}

int CashDisplay::get_n_col(const unsigned panel_ind) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||panel_ind < panels.size());
  }catch(GeneralError){
    std::cerr << "CashDisplay::get_n_col(): Error, out-of-bound panel_ind=" << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].n_col;
}

int CashDisplay::get_o_row(const unsigned panel_ind) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||panel_ind < panels.size());
  }catch(GeneralError){
    std::cerr << "CashDisplay::get_o_row(): Error, out-of-bound panel_ind=" << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].o_row;
}

int CashDisplay::get_o_col(const unsigned panel_ind) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||panel_ind < panels.size());
  }catch(GeneralError){
    std::cerr << "CashDisplay::get_o_row(): Error, out-of-bound panel_ind=" << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].o_col;
}

void CashDisplay::open_window(const std::string& window_title)
{
  char *cstr = new char[window_title.length() + 1];
  strcpy(cstr, window_title.c_str());
  OpenDisplay(cstr,window_row,window_col);
  window_is_open = true;
}

void CashDisplay::open_png(const std::string& directory_name)
{
  char *cstr = new char[directory_name.length() + 1];
  strcpy(cstr, directory_name.c_str());
  OpenPNG(cstr,window_row,window_col);
  png_is_open = true;
}

void CashDisplay::color_rgb(const unsigned char color_ind,const unsigned char r,const unsigned char g,const unsigned char b)
{
  if(window_is_open || png_is_open){
    std::cerr << "CashDisplay::color_rgb(): Error, call this method before openning window or png"  << std::endl;
    exit(-1);
  }

  ColorRGB(color_ind,r,g,b);
}

