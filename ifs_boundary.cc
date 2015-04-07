#include <algorithm>



//computes the angle and scale coordinates of this limit set.  
//
//this function finds the boundary of the limit set, expressed as 
//a bunch of words.  Then it splits the boundary into the f and g halves, 
//removes the prefix "f" from the f half, and finds where that interval is, 
//and how long it is.  The position gives theta and the length 
//gives lambda
bool ifs::compute_coordinates(double* theta, double* lambda, int n_depth) {
  
  int verbose = 0;
  
  //check if we are in the reasonable region 
  if (abs(z) > 1.0/sqrt(2.0) + 0.01) return false;
  
  //first, compute all the balls
  ifs temp_IFS;
  temp_IFS.set_params(z,z);
  temp_IFS.depth = n_depth;
  
  //find all the balls
  double min_r;
  if (!temp_IFS.minimal_enclosing_radius(min_r))  return false;
  if (!temp_IFS.circ_connected(min_r)) return false;
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  std::vector<Ball> balls(0);
  temp_IFS.compute_balls(balls, initial_ball, n_depth);
  
  //get a box which contains the balls
  cpx ll, ur;
  temp_IFS.box_containing_balls(balls, ll, ur);
  cpx box_center = 0.5*(ur + ll);
  double box_radius = 0.5*(ur.real() - ll.real());
  
  //make the box slightly larger so that we make sure the 
  //balls have some room around them
  box_radius *= 1.05;
  ll = cpx(box_center.real()-box_radius, box_center.imag()-box_radius);
  ur = cpx(box_center.real()+box_radius, box_center.imag()+box_radius);
  
  //figure out how big the pixels should be
  //(the ball radius ought to be about 2 pixel diameters)
  double desired_pixel_diameter = balls[0].radius/2.5;
  int num_pixels = int( (2*box_radius)/desired_pixel_diameter + 1 );
  if (num_pixels > 1000) num_pixels = 1000;  
  
  //create a trap grid from the balls
  TrapGrid TG;
  TG.reset_grid(ll, ur, num_pixels);
  TG.fill_pixels(balls);
  
  //compute the boundary (the third coordinate is unused, so 
  //we need to set it to zero)
  std::vector<Point3d<int> > pixel_boundary(0);
  TG.compute_boundary(pixel_boundary);
  for (int i=0; i<(int)pixel_boundary.size(); ++i) {
    pixel_boundary[i].z=0;
  }
  if (verbose>0) TG.show(NULL, &pixel_boundary, NULL, NULL, NULL);
  
  //get the list of uv words from the boundary
  std::vector<Bitword> unreduced_word_boundary(pixel_boundary.size());
  for (int i=0; i<(int)pixel_boundary.size(); ++i) {
    int ii = pixel_boundary[i].x;
    int jj = pixel_boundary[i].y;
    int ball_index = -1;
    if (TG.grid[ii][jj].z_ball_status > 0) {
      ball_index = TG.grid[ii][jj].closest_z_ball;
    } else {
      ball_index = TG.grid[ii][jj].closest_w_ball;
    }   
    unreduced_word_boundary[i] = Bitword( balls[ball_index].word, 
                                          balls[ball_index].word_len );
    //std::cout << i << ": " << unreduced_word_boundary[i] << "\n";
  }
  
  //reduce the boundary
  //this rotates the boundary so that it starts with a 0
  //and then removes any duplicates (which are next to each other)
  int start_index = 0;
  int M = (int)pixel_boundary.size();
  while (unreduced_word_boundary[start_index].reverse_get(0) != 1 || 
         unreduced_word_boundary[(start_index+1)%M].reverse_get(0) != 0) ++start_index;
  start_index = (start_index+1)%M;
  std::vector<Bitword> word_boundary(0);
  word_boundary.push_back(unreduced_word_boundary[start_index]);
  for (int i=1; i<M; ++i) {
    Bitword b = unreduced_word_boundary[(start_index+i)%M];
    if (b != word_boundary.back()) {
      word_boundary.push_back(b);
    }
  }
  
  if (verbose>0) { 
    for (int i=0; i<(int)word_boundary.size(); ++i) {
      std::cout << i << ": " << word_boundary[i] << "\n";
    }
  }
  
  //find where the boundary goes 0->1 or 1->0
  //note it'll have to switch 0->1 first
  //and there is a switch at the beginning of the list
  std::vector<int> s01(0);
  std::vector<int> s10(0);
  M = (int)word_boundary.size();
  s10.push_back(-1);
  for (int i=0; i<M; ++i) {
    int a = word_boundary[i].reverse_get(0);
    int b = word_boundary[(i+1)%M].reverse_get(0);
    if (a != b) {
      if (a == 0) {
        s01.push_back(i);
      } else {
        s10.push_back(i);
      }
    }
  }
  
  //find the largest block of 0's
  std::vector<int> gaps(s01.size());
  for (int i=0; i<(int)s01.size(); ++i) {
    gaps[i] = s01[i] - s10[i];
  }
  int largest_gap_ind = -1;
  for (int i=0; i<(int)s01.size(); ++i) {
    if (largest_gap_ind == -1 || gaps[largest_gap_ind] < gaps[i]) {
      largest_gap_ind = i;
    }
  }
  
  //get the largest block of 0's
  int zero_block_start = s10[largest_gap_ind]+1;
  int zero_block_end = s01[largest_gap_ind];
  int zero_block_middle = (zero_block_start + zero_block_end)/2;
  int zero_block_len = zero_block_end - zero_block_start+1;
  std::vector<Bitword> block0( word_boundary.begin() + zero_block_start,
                               word_boundary.begin() + zero_block_end + 1);
  
  //strip the 0's off
  std::vector<Bitword> stripped0(block0.size());
  for (int i=0; i<(int)block0.size(); ++i) {
    stripped0[i] = block0[i].suffix(block0[i].len-1);
  }
  
  if (verbose>0) {
    std::cout << "Zero block: " << zero_block_start << ", " 
                              << zero_block_middle << ", " 
                              << zero_block_end << ", length: " << zero_block_len << "\n";
  }
  
  //find where this stripped block exists in the word boundary
  int stripped_block_start=0;
  int stripped_block_end=0;
  int stripped_block_middle=0;
  int stripped_block_len=0;
  int wL = word_boundary[0].len;
  //find where the stripped interval begins and ends
  //to make sure we have the correct position, we 
  //need to find the one which is followed by its stripped follower
  //(before seeing another copy of the first one)
  std::vector<std::pair<Bitword, int> > just_starts(0);
  std::vector<std::pair<Bitword, int> > just_ends(0);
  Bitword a = stripped0[0];
  Bitword b = stripped0[1];
  Bitword y = stripped0[stripped0.size()-2];
  Bitword z = stripped0[stripped0.size()-1];
  for (int i=0; i<M; ++i) {
    Bitword t = word_boundary[i].prefix(wL-1);
    if (t == a || t == b) {
      just_starts.push_back( std::make_pair( t, i ) );
      //std::cout << "start: " << just_starts.back().first << " " << just_starts.back().second << "\n";
    }
    if (t == y || t == z) {
      just_ends.push_back( std::make_pair( t, i ) );
      //std::cout << "end: " << just_ends.back().first << " " << just_ends.back().second << "\n";
    }
  }
  for (int i=0; i<(int)just_starts.size(); ++i) {
    int ip1 = (i+1)%just_starts.size();
    if (just_starts[i].first == a && just_starts[ip1].first == b) {
      stripped_block_start = just_starts[i].second;
      break;
    }
  }
  for (int i=0; i<(int)just_ends.size(); ++i) {
    int ip1 = (i+1)%just_ends.size();
    if (just_ends[i].first == y && just_ends[ip1].first == z) {
      stripped_block_end = (just_ends[ip1].second)%M;
      break;
    }
  }
  
    
  if (stripped_block_end < stripped_block_start) {
    stripped_block_len = (stripped_block_end+M) - stripped_block_start + 1;
  } else if (stripped_block_end > stripped_block_start) {
    stripped_block_len = stripped_block_end - stripped_block_start + 1;
  } else if (stripped_block_end == stripped_block_start) {
    stripped_block_len = M;
  }
  stripped_block_middle = (stripped_block_start + stripped_block_len/2)%M;
  
  if (verbose>0) {
    std::cout << "Stripped block: " << stripped_block_start << ", " 
                                    << stripped_block_middle << ", " 
                                  << stripped_block_end << ", length: " << stripped_block_len << "\n";
  }
  
  if (stripped_block_middle >= zero_block_middle) {
    *theta = double(stripped_block_middle - zero_block_middle)/double(M);
  } else {
    *theta = double( (stripped_block_middle+M) - zero_block_middle)/double(M);
  }
  if (*theta > 0.5) *theta -= 1.0;
  if (verbose>0) std::cout << "theta: " << *theta << "\n";
  
  *lambda = double(stripped_block_len)/double(zero_block_len);
  if (verbose>0) std::cout << "lambda: " << *lambda << "\n";
  
  
  return true;
}






//this computes only theta with the new method
bool ifs::compute_new_theta(double* theta, int n_depth) {
  
  int verbose = 0;
  
  std::vector<Bitword> word_boundary;
  std::vector<Bitword> zero_word_boundary;
  
  if (!compute_boundary_and_f_boundary(word_boundary, 
                                       zero_word_boundary, 
                                       n_depth,
                                       verbose)) {
    return false;
  }
  
  //*****************************************************************
  
  //find where the first two and last two zero words in the big list are in the 
  //zero list
  Bitword z0 = word_boundary[0];
  Bitword z1 = word_boundary[1];
  int last_zero=0;
  while (word_boundary[last_zero+1].reverse_get(0) == 0) ++last_zero;
  if (last_zero==0) {
    std::cout << "Weird parameter: " << z << "\n";
    return false;
  }
  Bitword z_1 = word_boundary[last_zero];
  Bitword z_2 = word_boundary[last_zero-1];
 
  if (verbose>0) {
    std::cout << "Found first two zero words: " << z0 << " " << z1 << "\n";
    std::cout << "And last two: " << z_2 << " " << z_1 << "\n";
  }
  
  //figure out where the block of main zeros is in the list of zeros
  
  int zero_block_start=0;
  int zero_block_last=0;
  int M = (int)zero_word_boundary.size();
  for (int i=0; i<(int)zero_word_boundary.size(); ++i) {
    if (zero_word_boundary[i] == z0 && zero_word_boundary[(i+1)%M] == z1) {
      zero_block_start = i;
    }
    if (zero_word_boundary[i] == z_2 && zero_word_boundary[(i+1)%M] == z_1) {
      zero_block_last = (i+1)%M;
    }
  }
  
  if (verbose > 0) {
    std::cout << "Zero block start: " << zero_block_start << "\n";
    std::cout << "Zero block last: " << zero_block_last << "\n";
  }
  int zero_block_len = -1;
  if (zero_block_start == zero_block_last || 
      (zero_block_last+1)%M == zero_block_start) {
    zero_block_len = zero_word_boundary.size();
  } else {
    zero_block_len = (zero_block_last+1)-zero_block_start;
  }
  if (zero_block_len < 0) zero_block_len += M;
  int zero_block_remainder = M - zero_block_len;
  
  if (verbose > 0) {
    std::cout << "zero block len: " << zero_block_len << "\n";
    std::cout << "Zero block remainder: " << zero_block_remainder << "\n";
  }
  
  //find the position of the stripped middle zero block entry
  int zero_block_middle = (zero_block_start + zero_block_len/2)%M;
  int zero_block_middle_in_word_boundary = 0 + zero_block_len/2;
  Bitword stripped_middle = zero_word_boundary[zero_block_middle];
  stripped_middle = stripped_middle.suffix(stripped_middle.len-1);
  
  int stripped_middle_preimage = 0;
  for (int i=0; i<(int)word_boundary.size(); ++i) {
    if (word_boundary[i].prefix(word_boundary[i].len-1) == stripped_middle) {
      stripped_middle_preimage = i;
      break;
    }
  }
  
  if (verbose>0) {
    std::cout << "zero block middle: " << zero_block_middle << "\n";
    std::cout << "zero block middle in word boundary: " << zero_block_middle_in_word_boundary << "\n";
    std::cout << "Stripped middle word: " << stripped_middle << "\n";
    std::cout << "Stripped middle preimage: " << stripped_middle_preimage << "\n";
  }
  
  double distance_to_preimage = 0;
  double distance_completely_around = 0;
  double maximal_insertion_amount = zero_block_remainder;
  bool passed_preimage = false;
  M = word_boundary.size();
  int raw_distance_to_preimage = 0;
  int raw_distance_around = 0;
  for (int i=(zero_block_middle_in_word_boundary+1)%M; 
           i!=zero_block_middle_in_word_boundary; 
           i=(i+1)%M) {
    if (i==stripped_middle_preimage) passed_preimage = true;
    int common_prefix = word_boundary[i].common_prefix(word_boundary[(i+1)%M]);
    double amount_from_this_step = pow(abs(z), common_prefix)*maximal_insertion_amount;
    if (!passed_preimage) {
      ++raw_distance_to_preimage;
      distance_to_preimage += amount_from_this_step+1;
    }
    ++raw_distance_around;
    distance_completely_around += amount_from_this_step+1;
  }
  
  if (verbose>0) {
    std::cout << "Maximal insertion amount: " << maximal_insertion_amount << "\n";
    std::cout << "Raw distance to preimage: " << raw_distance_to_preimage << "\n";
    std::cout << "Distance to preimage: " << distance_to_preimage << "\n";
    std::cout << "Raw distance around: " << raw_distance_around << "\n";
    std::cout << "Distance completely around: " << distance_completely_around << "\n";
  }
  
  *theta = distance_to_preimage / distance_completely_around;
  return true;
}
  





//returns the complete boundary of the limit set, in words
//and the complete boundary of fL.  They are shifted so they begin
//with as many zeros as possible (i.e. fL begins 00)
bool ifs::compute_boundary_and_f_boundary(std::vector<Bitword>& whole_boundary, 
                                          std::vector<Bitword>& f_boundary,
                                          int n_depth,
                                          int verbose) {
  
  //check if we are in the reasonable region 
  if (abs(z) > 1.0/sqrt(2.0) + 0.01) return false;
  
  //first, compute all the balls
  ifs temp_IFS;
  temp_IFS.set_params(z,z);
  temp_IFS.depth = n_depth;
  
  //find all the balls
  double min_r;
  if (!temp_IFS.minimal_enclosing_radius(min_r))  return false;
  if (!temp_IFS.circ_connected(min_r)) return false;
  
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  std::vector<Ball> balls(0);
  temp_IFS.compute_balls(balls, initial_ball, n_depth);
  
  //get a box which contains the balls
  cpx ll, ur;
  temp_IFS.box_containing_balls(balls, ll, ur);
  cpx box_center = 0.5*(ur + ll);
  double box_radius = 0.5*(ur.real() - ll.real());
  
  //make the box slightly larger so that we make sure the 
  //balls have some room around them
  box_radius *= 1.05;
  ll = cpx(box_center.real()-box_radius, box_center.imag()-box_radius);
  ur = cpx(box_center.real()+box_radius, box_center.imag()+box_radius);
  
  //figure out how big the pixels should be
  //(the ball radius ought to be about 2 pixel diameters)
  double desired_pixel_diameter = balls[0].radius/2.5;
  int num_pixels = int( (2*box_radius)/desired_pixel_diameter + 1 );
  if (num_pixels > 1000) num_pixels = 1000;  
  
  //create a trap grid from the balls
  TrapGrid TG;
  TG.reset_grid(ll, ur, num_pixels);
  TG.fill_pixels(balls);
  
  //compute the boundary (the third coordinate is unused, so 
  //we need to set it to zero)
  std::vector<Point3d<int> > pixel_boundary(0);
  TG.compute_boundary(pixel_boundary);
  for (int i=0; i<(int)pixel_boundary.size(); ++i) {
    pixel_boundary[i].z=0;
  }
  if (verbose>0) TG.show(NULL, &pixel_boundary, NULL, NULL, NULL);
  
  //get the list of uv words from the boundary
  std::vector<Bitword> unreduced_word_boundary(pixel_boundary.size());
  for (int i=0; i<(int)pixel_boundary.size(); ++i) {
    int ii = pixel_boundary[i].x;
    int jj = pixel_boundary[i].y;
    int ball_index = -1;
    if (TG.grid[ii][jj].z_ball_status > 0) {
      ball_index = TG.grid[ii][jj].closest_z_ball;
    } else {
      ball_index = TG.grid[ii][jj].closest_w_ball;
    }   
    unreduced_word_boundary[i] = Bitword( balls[ball_index].word, 
                                          balls[ball_index].word_len );
    //std::cout << i << ": " << unreduced_word_boundary[i] << "\n";
  }
  
  if (verbose>0) {
    std::cout << "Unreduced word boundary:\n";
    for (int i=0; i<(int)unreduced_word_boundary.size(); ++i) {
      std::cout << i << ": " << unreduced_word_boundary[i] << "\n";
    }
    std::cout << "\n";
  }
  
  //reduce the boundary
  //this rotates the boundary so that it starts with a 0
  //and then removes any duplicates (which are next to each other)
  int start_index = 0;
  int M = (int)pixel_boundary.size();
  while (unreduced_word_boundary[start_index].reverse_get(0) != 1 || 
         unreduced_word_boundary[(start_index+1)%M].reverse_get(0) != 0) ++start_index;
  start_index = (start_index+1)%M;
  std::vector<Bitword> word_boundary(0);
  word_boundary.push_back(unreduced_word_boundary[start_index]);
  for (int i=1; i<M; ++i) {
    Bitword b = unreduced_word_boundary[(start_index+i)%M];
    if (b != word_boundary.back()) {
      word_boundary.push_back(b);
    }
  }
  
  if (verbose>0) { 
    std::cout << "Full word boundary:\n";
    for (int i=0; i<(int)word_boundary.size(); ++i) {
      std::cout << i << ": " << word_boundary[i] << "\n";
    }
  }



  //**********************************************************

  //now we have to repeat this with just the 0's
  //we'll use the same size grid to make sure 
  //all the balls match up
  
  //get the zero balls
  std::vector<Ball> zero_balls(0);
  for (int i=0; i<(int)balls.size(); ++i) {
    Bitword b (balls[i].word, balls[i].word_len);
    if (b.reverse_get(0) == 0) {
      zero_balls.push_back(balls[i]);
    }
  }
  
  //create a trap grid from the balls
  TrapGrid TG0;
  TG0.reset_grid(ll, ur, num_pixels);
  TG0.fill_pixels(zero_balls);
  
  //compute the boundary (the third coordinate is unused, so 
  //we need to set it to zero)
  std::vector<Point3d<int> > zero_pixel_boundary(0);
  TG0.compute_boundary(zero_pixel_boundary);
  for (int i=0; i<(int)zero_pixel_boundary.size(); ++i) {
    pixel_boundary[i].z=0;
  }
  if (verbose>0) TG0.show(NULL, &zero_pixel_boundary, NULL, NULL, NULL);

  //get the list of uv words from the boundary
  std::vector<Bitword> zero_unreduced_word_boundary(zero_pixel_boundary.size());
  for (int i=0; i<(int)zero_pixel_boundary.size(); ++i) {
    int ii = zero_pixel_boundary[i].x;
    int jj = zero_pixel_boundary[i].y;
    int ball_index = -1;
    if (TG0.grid[ii][jj].z_ball_status > 0) {
      ball_index = TG0.grid[ii][jj].closest_z_ball;
    } else {
      ball_index = TG0.grid[ii][jj].closest_w_ball;
    }   
    zero_unreduced_word_boundary[i] = Bitword( zero_balls[ball_index].word, 
                                               zero_balls[ball_index].word_len );
    //std::cout << i << ": " << unreduced_word_boundary[i] << "\n";
  }
  
  //remove any duplicates (which are next to each other)
  M = (int)zero_unreduced_word_boundary.size();
  std::vector<Bitword> zero_word_boundary(0);
  zero_word_boundary.push_back(zero_unreduced_word_boundary[0]);
  for (int i=1; i<M; ++i) {
    Bitword b = zero_unreduced_word_boundary[i];
    if (b != zero_word_boundary.back()) {
      zero_word_boundary.push_back(b);
    }
  }
  while (zero_word_boundary.back() == zero_word_boundary.front()) {
    zero_word_boundary.erase( zero_word_boundary.begin() + (zero_word_boundary.size()-1) );
  }
  
  //rotate so it starts 00 and ends 01
  int k=0;
  M = (int)zero_word_boundary.size();
  while (zero_word_boundary[k].reverse_get(1) == 0) k = (k+1)%M;
  while (zero_word_boundary[k].reverse_get(1) == 1) k = (k+1)%M;
  //k is now the index of the first place we have 00
  std::vector<Bitword> new_vec(zero_word_boundary.begin() + k, zero_word_boundary.end());
  new_vec.insert(new_vec.end(), zero_word_boundary.begin(), zero_word_boundary.begin()+k);
  zero_word_boundary = new_vec;
  
  
  if (verbose>0) { 
    std::cout << "Zero word boundary:\n";
    for (int i=0; i<(int)zero_word_boundary.size(); ++i) {
      std::cout << i << ": " << zero_word_boundary[i] << "\n";
    }
  }
  
  whole_boundary = word_boundary;
  f_boundary = zero_word_boundary;
  
  return true;
}



bool ifs::compute_boundary_space(BoundarySpace& BS,
                                 int n_depth,
                                 int lam_depth) {
  int verbose = 0;
  
  std::vector<Bitword> word_boundary;
  std::vector<Bitword> zero_word_boundary;
  
  if (!compute_boundary_and_f_boundary(word_boundary, 
                                       zero_word_boundary, 
                                       n_depth,
                                       verbose)) {
    return false;
  }
  
  if (verbose>0) {
    std::cout << "Whole boundary: \n";
    for (int i=0; i<(int)word_boundary.size(); ++i) {
      std::cout << i << ": " << word_boundary[i] << "\n";
    }
    std::cout << "Zero boundary: \n";
    for (int i=0; i<(int)zero_word_boundary.size(); ++i) {
      std::cout << i << ": " << zero_word_boundary[i] << "\n";
    }
  }
  
  //find the interval in the zero boundary which lies on the 
  //inside of the circle
  std::vector<Bitword> interior_zeros(0);
  Bitword initial_zeros_1 = word_boundary[0];
  Bitword initial_zeros_2 = word_boundary[1];
  int ones_position = 0;
  while (word_boundary[ones_position].reverse_get(0) == 0) ++ones_position;
  Bitword final_zeros_1 = word_boundary[ones_position-2];
  Bitword final_zeros_2 = word_boundary[ones_position-1];
  int M = (int)zero_word_boundary.size();
  for (int i=0; i<M; ++i) {
    if (zero_word_boundary[i] == final_zeros_1 &&
        zero_word_boundary[(i+1)%M] == final_zeros_2) {
      for (int j=2; j<M; ++j) {
        if (zero_word_boundary[(i+j)%M] == initial_zeros_1 &&
            zero_word_boundary[(i+j+1)%M] == initial_zeros_2) break;
        interior_zeros.push_back(zero_word_boundary[(i+j)%M]);
      }
      break;
    }
  }
  
  if (verbose>0) {
    std::cout << "Ones position: " << ones_position << "\n";
    std::cout << "Initial zeros: " << initial_zeros_1 << " " << initial_zeros_2 << "\n";
    std::cout << "Final zeros: " << final_zeros_1 << " " << final_zeros_2 << "\n";
    std::cout << "Interior zeros:\n";
    for (int i=0; i<(int)interior_zeros.size(); ++i) {
      std::cout << i << ": " << interior_zeros[i] << "\n";
    }
  }
  
  //repeat the same thing for the ones
  std::vector<Bitword> interior_ones(0);
  Bitword initial_ones_1 = word_boundary[ones_position];
  Bitword initial_ones_2 = word_boundary[(ones_position+1)%word_boundary.size()];
  Bitword final_ones_1 = word_boundary[word_boundary.size()-2];
  Bitword final_ones_2 = word_boundary[word_boundary.size()-1];
  //swap the first bit so it's easy
  initial_ones_1.reverse_set(0,0); initial_ones_2.reverse_set(0,0);
  final_ones_1.reverse_set(0,0); final_ones_2.reverse_set(0,0);
  for (int i=0; i<M; ++i) {
    if (zero_word_boundary[i] == final_ones_1 &&
        zero_word_boundary[(i+1)%M] == final_ones_2) {
      for (int j=2; j<M; ++j) {
        if (zero_word_boundary[(i+j)%M] == initial_ones_1 &&
            zero_word_boundary[(i+j+1)%M] == initial_ones_2) break;
        interior_ones.push_back(zero_word_boundary[(i+j)%M]);
        interior_ones.back().reverse_set(0,1);
      }
      break;
    }
  }
  
  if (interior_ones.size() == 0) { //we were unlucky and didn't find the 1's boundary
    return false;
  }
  
  if (verbose>0) {
    std::cout << "Initial ones (first digit swapped): " << initial_ones_1 << " " << initial_ones_2 << "\n";
    std::cout << "Final ones: " << final_ones_1 << " " << final_ones_2 << "\n";
    std::cout << "Interior ones:\n";
    for (int i=0; i<(int)interior_ones.size(); ++i) {
      std::cout << i << ": " << interior_ones[i] << "\n";
    }
  }
  
  
  std::vector<std::vector<Bitword> > level_interior_zeros(lam_depth);
  std::vector<std::vector<Bitword> > level_interior_ones(lam_depth);
  level_interior_zeros[0] = interior_zeros;
  level_interior_ones[0] = interior_ones;
  for (int i=1; i<lam_depth; ++i) {
    level_interior_zeros[i].resize(1);
    level_interior_zeros[i][0] = interior_zeros[0].shift_right(i, 0);
    for (int j=1; j<(int)interior_zeros.size(); ++j) {
      Bitword next = interior_zeros[j].shift_right(i, 0);
      if (next != level_interior_zeros[i].back()) {
        level_interior_zeros[i].push_back(next);
      }
    }
    
    level_interior_ones[i].resize(1);
    level_interior_ones[i][0] = interior_ones[0].shift_right(i, 1);
    for (int j=1; j<(int)interior_ones.size(); ++j) {
      Bitword next = interior_ones[j].shift_right(i, 1);
      if (next != level_interior_ones[i].back()) {
        level_interior_ones[i].push_back(next);
      }
    }
  }
  
  if (verbose>0) {
    for (int i=1; i<lam_depth; ++i) {
      std::cout << "level " << i << " interior zeros:\n";
      for (int j=0; j<(int)level_interior_zeros[i].size(); ++j) {
        std::cout << j << ": " << level_interior_zeros[i][j] << "\n";
      }
      std::cout << "level " << i << " interior ones:\n";
      for (int j=0; j<(int)level_interior_ones[i].size(); ++j) {
        std::cout << j << ": " << level_interior_ones[i][j] << "\n";
      }
    }
  }
  
  //create the boundary
  BS.cut_depth = n_depth;
  BS.lam_depth = lam_depth;
  BS.boundary = word_boundary;
  BS.b_word_depth = std::vector<int>(word_boundary.size(), 0);
  BS.lam.resize(0);
  
  std::set<Bitword> boundary_set(word_boundary.begin(), word_boundary.end());
  
  for (int cut_level=0; cut_level < lam_depth; ++cut_level) {
    
    if (verbose>0) {
      std::cout << "Boundary before cutting at level " << cut_level << ":\n";
      for (int i=0; i<(int)BS.boundary.size(); ++i) {
        std::cout << i << ": " << BS.boundary[i] << "\n";
      }
      std::cout << "Lamination: ";
      for (int i=0; i<(int)BS.lam.size(); ++i) {
        std::cout << "(" << BS.lam[i].x << "," << BS.lam[i].y << "," << BS.lam[i].z << ")";
      }
      std::cout << "\n";
    }    
    
    //find all the places in the boundary where 
    //the highest swap place is at level cut_level
    M = (int)BS.boundary.size();
    std::vector<int> swaps(0);
    if (BS.boundary.back().common_prefix(BS.boundary[0]) == cut_level) {
      swaps.push_back(0);
    }
    for (int i=0; i<(int)BS.boundary.size()-1; ++i) {
      if (BS.boundary[i].common_prefix(BS.boundary[i+1]) == cut_level) {
        swaps.push_back(i+1);
      }
    }
    
    if (verbose>0) {
      std::cout << "Swaps: ";
      for (int i=0; i<(int)swaps.size(); ++i) {
        std::cout << swaps[i] << " ";
      }
      std::cout << "\n";
    }
    
    //go *backward* through the swaps so that we don't screw up the indices
    for (int i=(int)swaps.size()-1; i>=0; --i) {
      std::vector<Bitword> list_to_insert;
      if (BS.boundary[swaps[i]].reverse_get(cut_level) == 0) {
        list_to_insert = level_interior_zeros[cut_level];
      } else {
        list_to_insert = level_interior_ones[cut_level];
      }
      for (int j=0; j<(int)list_to_insert.size(); ++j) {
        list_to_insert[j].copy_prefix(BS.boundary[swaps[i]], cut_level);
      }
      std::vector<Bitword> trimmed_list_to_insert(0);
      for (int j=0; j<(int)list_to_insert.size(); ++j) {
        if ( boundary_set.find(list_to_insert[j]) == boundary_set.end() ) {
          trimmed_list_to_insert.push_back(list_to_insert[j]);
          boundary_set.insert(list_to_insert[j]);
        }
      }
      int i_size = trimmed_list_to_insert.size();
      
      BS.boundary.insert(BS.boundary.begin() + swaps[i], 
                         trimmed_list_to_insert.begin(), 
                         trimmed_list_to_insert.end());
      BS.b_word_depth.insert(BS.b_word_depth.begin() + swaps[i],
                             i_size, cut_level+1);
      for (int j=0; j<(int)BS.lam.size(); ++j) {
        if (BS.lam[j].x >= swaps[i]) BS.lam[j].x += i_size;
        if (BS.lam[j].y >= swaps[i]) BS.lam[j].y += i_size;
      }
      //note that swaps[i] now points to where the list starts
      //let's keep the higher swaps correct
      for (int j=i+1; j<(int)swaps.size(); ++j) {
        swaps[j] += i_size;
      }
    }
    
    //now add in the new leaves; we add a leaf whenever they have the same 
    //prefix of *length* cut_level
    for (int i=0; i<(int)swaps.size(); ++i) {
      for (int j=i+1; j<(int)swaps.size(); ++j) {
        if (BS.boundary[swaps[i]].common_prefix(BS.boundary[swaps[j]]) == cut_level) {
          BS.lam.push_back(Point3d<int>(swaps[i], swaps[j], cut_level));
        }
      }
    }
      
  }
  
  return true;
}





/***************************************************************************
 * these are plain ball functions which are useful for building
 * the non-grid ball boundary
 * *************************************************************************/

/****************************************************************************
 * get a list of all the indices of balls which intersect the ball b.
 * These indices are returned as a list of ints, where each int gives 
 * the fg word of the intersecting ball in bits
 * **************************************************************************/
std::vector<int> balls_which_intersect_ball(Ball& b, 
                                            std::vector<std::pair<Ball,Ball> >& intersection_pairs,
                                            int verbose) {
  std::vector<Bitword> bitword_ans(0);
  std::vector<int> ans;
  //this is the same as the uv graph
  //we just need to find all suffixes of b which agree with a prefix of an 
  //intersection pair and then swap for the other member of the 
  //intersection pair
  Bitword bw(b.word, b.word_len);
  
  //get a list of the bitwords from the balls
  std::vector<std::pair<Bitword,Bitword> > intersection_pair_bitwords(intersection_pairs.size());
  for (int i=0; i<(int)intersection_pairs.size(); ++i) {
    intersection_pair_bitwords[i] = std::make_pair(Bitword(intersection_pairs[i].first.word, intersection_pairs[i].first.word_len),
                                                   Bitword(intersection_pairs[i].second.word, intersection_pairs[i].second.word_len));
  }
  
  //do the swapping
  for (int n=1; n<=bw.len; ++n) {
    Bitword suffix = bw.suffix(n);
    for (int i=0; i<(int)intersection_pairs.size(); ++i) {
      if (suffix.reverse_get(0) == 0) {
        if (intersection_pair_bitwords[i].first == suffix) {
          bitword_ans.push_back( bw.swapped_suffix(n, intersection_pair_bitwords[i].second) );
        }
      } else {
        if (intersection_pair_bitwords[i].second == suffix) {
          bitword_ans.push_back( bw.swapped_suffix(n, intersection_pair_bitwords[i].first) );
        }
      }
    }
  }
  
  ans.resize(bitword_ans.size());
  for (int i=0; i<(int)bitword_ans.size(); ++i) {
    ans[i] = bitword_ans[i].to_int();
  }
  
  return ans;
}


/***************************************************************************
 * this returns true if the three given angles are cyclically in order
 * *************************************************************************/
bool cyclically_ordered(double a, double b, double c) {
  double PI = 3.14159265358979;
  if (b<a) b += 2*PI;
  if (c<a) c += 2*PI;
  return b<c;
}

/***************************************************************************
 * this searches among the indices in intersecting balls to find the 
 * one which intersects last before the angle.  If the angle is inside 
 * an intersecting interval, the result is undefined
 * *************************************************************************/
int last_intersecting_ball_before_angle(std::vector<Ball>& balls, 
                                        int current_ball_i, 
                                        std::vector<int>& intersecting_balls, 
                                        double angle,
                                        int verbose) {
  if (verbose>0) {
    std::cout << "Looking amongst intersecting balls for the first one after angle " << angle << "\n";
    for (int i=0; i<(int)intersecting_balls.size(); ++i) {
      std::cout << intersecting_balls[i] << " ";
    } std::cout << "\n";
  }
  std::vector<double> last_angles(intersecting_balls.size());
  for (int i=0; i<(int)intersecting_balls.size(); ++i) {
    last_angles[i] = balls[current_ball_i].intersection_interval( balls[intersecting_balls[i]] ).second;
  }
  int closest_last_angle_i = 0;
  for (int i=1; i<(int)last_angles.size(); ++i) {
    if (cyclically_ordered( last_angles[closest_last_angle_i], last_angles[i], angle )) {
      closest_last_angle_i = i;
    }
  }
    
  return intersecting_balls[closest_last_angle_i];
}

/**************************************************************************
 * this searchs among the indices in intersecting_balls to find the 
 * index of the ball which intersects current_ball_i first after the 
 * ball prev_ball_i
 * ************************************************************************/
int next_intersecting_ball_after_ball(std::vector<Ball>& balls, 
                                      int current_ball_i, 
                                      int prev_ball_i, 
                                      std::vector<int>& intersecting_balls,
                                      int verbose) {
  if (verbose>0) {
    std::cout << "Looking amongst intersecting balls for the first one after " << prev_ball_i << "\n";
    for (int i=0; i<(int)intersecting_balls.size(); ++i) {
      std::cout << intersecting_balls[i] << " ";
    } std::cout << "\n";
  }
  double prev_angle = balls[current_ball_i].intersection_interval( balls[prev_ball_i] ).second;
  std::vector<double> first_angles(intersecting_balls.size());
  for (int i=0; i<(int)intersecting_balls.size(); ++i) {
    first_angles[i] = balls[current_ball_i].intersection_interval( balls[intersecting_balls[i]] ).first;
  }
  int closest_first_angle_i = 0;
  for (int i=1; i<(int)first_angles.size(); ++i) {
    if (cyclically_ordered( prev_angle, first_angles[i], first_angles[closest_first_angle_i] )) {
      closest_first_angle_i = i;
    }
  }
  return intersecting_balls[closest_first_angle_i];
}







/****************************************************************************
* compute boundary balls not using a trap grid
* this function will always start with the first (top) f-ball
* *unless* only_1_section is true, in which case it'll return 
* the set of balls on the boundary from the rightmost 1-ball to the junction 
* with the 0 side
****************************************************************************/
void non_grid_ball_boundary_indices(std::vector<int>& boundary, 
                                    std::vector<Ball>& balls, 
                                    std::vector<std::pair<Ball,Ball> >& intersection_pairs,
                                    bool only_1_section,
                                    int verbose) {
  boundary.resize(0);
  int ball_depth = balls[0].word_len;
  
  //find the ball with the greatest x coord 
  int max_x_ball_i = 0;
  for (int i=0; i<(int)balls.size(); ++i) {
    if (balls[i].center.real() > balls[max_x_ball_i].center.real()) {
      max_x_ball_i = i;
    }
  }
  
  if (verbose>0) {
    std::cout << "Finding boundary ball indices\n";
    std::cout << "Starting ball (" << max_x_ball_i << "," << Bitword(max_x_ball_i, ball_depth) << "): " << balls[max_x_ball_i] << "\n";
  } 
  
  //so balls[max_x_ball_i] must have angle 0 exposed
  int current_ball_i = max_x_ball_i;
  std::vector<int> intersecting_balls = balls_which_intersect_ball(balls[current_ball_i], intersection_pairs);
  //get the ball which is previous in the boundary
  int prev_ball_i = last_intersecting_ball_before_angle(balls, current_ball_i, intersecting_balls, 0.0,verbose);  
  int next_ball_i;
  
  while (true) {
  
    if (verbose>0) {
      std::cout << "Prev ball: " << prev_ball_i << ", " << Bitword(prev_ball_i, ball_depth) << "\n";
      std::cout << "Current ball: " << current_ball_i << ", " << Bitword(current_ball_i, ball_depth) << "\n";
    }
    
    //get the next ball after the previous one
    next_ball_i = next_intersecting_ball_after_ball(balls, current_ball_i, prev_ball_i, intersecting_balls, verbose);
    
    if (verbose>0) {
      std::cout << "Next ball: " << next_ball_i << ", " << Bitword(next_ball_i, ball_depth) << "\n";
    }
    
    //we might have looped around
    if (boundary.size() > 2 &&
        current_ball_i == boundary[0] &&
        next_ball_i == boundary[1]) {
      if (verbose>0) {
        std::cout << "We've looped\n";
      }
      break;
    }
    
    //if we haven't looped, we should put on the current ball
    boundary.push_back(current_ball_i);
    
    //next becomes current and current becomes previous
    prev_ball_i = current_ball_i;
    current_ball_i = next_ball_i;
    
    //if we're suppose to bail out when we hit the 0 side, we need to check 
    //for that
    if (only_1_section) {
      if ( ((current_ball_i >> (ball_depth-1))&1) == 0 ) return;
    }
    
    //compute the intersecting balls
    intersecting_balls = balls_which_intersect_ball(balls[current_ball_i], intersection_pairs, verbose);
  }
  
  //we *always* start with a g-ball (maximal x-coord), 
  //so it suffices to scan until we find an f-ball, and rotate the 
  //list to start there
  int f_ball_i=0;
  while (balls[boundary[f_ball_i]].last_gen_index() == 1) f_ball_i++;
  
  std::rotate(boundary.begin(), boundary.begin() + f_ball_i, boundary.end());
  
}


/***************************************************************************
 * similar to the above function, except that it only finds the pair 
 * of balls which transition from the g side to the f side (on top)
 * *************************************************************************/
void non_grid_ball_boundary_transition_indices(Bitword& g_bitword,
                                               Bitword& f_bitword,
                                               std::vector<Ball>& balls, 
                                               std::vector<std::pair<Ball,Ball> >& intersection_pairs,
                                               int verbose) {
  int ball_depth = balls[0].word_len;
  
  //find the ball with the greatest x coord 
  int max_x_ball_i = 0;
  for (int i=0; i<(int)balls.size(); ++i) {
    if (balls[i].center.real() > balls[max_x_ball_i].center.real()) {
      max_x_ball_i = i;
    }
  }
  
  if (verbose>0) {
    std::cout << "Finding boundary ball indices\n";
    std::cout << "Starting ball (" << max_x_ball_i << "," << Bitword(max_x_ball_i, ball_depth) << "): " << balls[max_x_ball_i] << "\n";
  } 
  
  //so balls[max_x_ball_i] must have angle 0 exposed
  int current_ball_i = max_x_ball_i;
  std::vector<int> intersecting_balls = balls_which_intersect_ball(balls[current_ball_i], intersection_pairs);
  //get the ball which is previous in the boundary
  int prev_ball_i = last_intersecting_ball_before_angle(balls, current_ball_i, intersecting_balls, 0.0,verbose);  
  int next_ball_i;
  
  while (true) {
  
    if (verbose>0) {
      std::cout << "Prev ball: " << prev_ball_i << ", " << Bitword(prev_ball_i, ball_depth) << "\n";
      std::cout << "Current ball: " << current_ball_i << ", " << Bitword(current_ball_i, ball_depth) << "\n";
    }
    
    //get the next ball after the previous one
    next_ball_i = next_intersecting_ball_after_ball(balls, current_ball_i, prev_ball_i, intersecting_balls, verbose);
    
    if (verbose>0) {
      std::cout << "Next ball: " << next_ball_i << ", " << Bitword(next_ball_i, ball_depth) << "\n";
    }
    
    //we should put on the current ball
    //** we don't have to do this because we're just getting the pair**
    //boundary.push_back(current_ball_i);
    
    //next becomes current and current becomes previous
    prev_ball_i = current_ball_i;
    current_ball_i = next_ball_i;
    
    //if the previous ball starts with 1 and the current ball 
    //starts with 0, then we need to just return these two numbers
    if ( ((current_ball_i >> (ball_depth-1))&1) == 0 &&
         ((prev_ball_i >> (ball_depth-1))&1) == 1) {
      g_bitword = Bitword(prev_ball_i, ball_depth);
      f_bitword = Bitword(current_ball_i, ball_depth);
      return;
    }
    
    //compute the intersecting balls
    intersecting_balls = balls_which_intersect_ball(balls[current_ball_i], intersection_pairs, verbose);
  }
  
}









/**************************************************************************
 * find the peak of a discrete function
 * ************************************************************************/
void find_integer_peak(std::vector<std::pair<int,bool> >& L, int& p1, int& p2, int& p3) {
  //pare the list down
  std::vector<std::pair<int,int> > true_L(0);
  for (int i=0; i<(int)L.size(); ++i) {
    if (L[i].second) {
      true_L.push_back(std::make_pair(L[i].first, i));
    }
  }
  int last_increase = -1;
  int last_decrease = -1;
  for (int i=0; i<(int)true_L.size()-1; ++i) {
    if ( true_L[i].first > true_L[i+1].first ) {
      last_decrease = i;
      if (last_increase >= 0) {
        p1 = true_L[last_increase].second;
        //std::cout << last_increase << " " << last_decrease << " " << last_increase + 1 + (last_decrease - last_increase)/2;
        p2 = true_L[ last_increase + 1 + (last_decrease - last_increase)/2 ].second;
        p3 = true_L[last_decrease+1].second;
        return;
      }
    } else if ( true_L[i].first < true_L[i+1].first ) {
      last_increase = i;
    }
  }
  p1 = p2 = p3 = -1;
  return;
}

/*************************************************************************
 * find the closest ball which doesn't share a prefix with the specified
 * ball of the given length
 *************************************************************************/
double ifs::matching_prefix_radius(Bitword& b, 
                                   int prefix_len) {
  //suppose that the ball starts with 1
  //then we need to find the distance to any 0 ball
  //then recursively we need to find the distance 
  //between this ball shifted left (wlog a 0 ball) and the 1 balls, etc
  //std::cout << "Finding distance from " << b << " and anything without prefix length " << prefix_len << "\n";
  if (prefix_len == 1) {
    return distance_from_other_half(b);
  } else {
    double d1 = distance_from_other_half(b);
    Bitword b2 = b.suffix(b.len-1);
    double d2 = abs(z)*matching_prefix_radius(b2, prefix_len-1);
    double d = (d1 < d2 ? d1 : d2);
    return d;
  }
}

/*************************************************************************
 * get the distance from this ball to any ball with a different first letter
 *************************************************************************/
double ifs::distance_from_other_half(Bitword& bw) {
  double min_r;
  minimal_enclosing_radius(min_r);
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  Ball test_ball = apply_bitword(bw, initial_ball);
  std::deque<Ball> stack(0);
  stack.push_back( act_on_left( 1-bw.reverse_get(0), initial_ball ) );
  std::vector<double> min_dists(0);
  //get the initial distance
  //std::cout << "Getting distance to other side of " << bw << "\n";
  min_dists.push_back( (abs(stack.back().center - test_ball.center) - stack.back().radius) - test_ball.radius );
  while (stack.size() > 0) {
    Ball b = stack.front();
    stack.pop_front();
    //get the distance
    double d = (abs(b.center-test_ball.center) - b.radius) - test_ball.radius;
    if ( (int)min_dists.size() >= b.word_len ) { //there is some distance for a ball of this depth
      if (b.word_len < test_ball.word_len &&  d < min_dists[b.word_len-1] + 2*b.radius ) {
        stack.push_back( act_on_right(0, b) );
        stack.push_back( act_on_right(1, b) );
      }
      if ( d < min_dists[b.word_len-1] ) {
        min_dists[b.word_len-1] = d;
        //std::cout << "New min dist of " << d << " with " << Bitword(b.word, b.word_len) << "\n";
      }
    } else { //there isn't a distance for a ball of this depth, so we need to put a new one
      if ( d < min_dists[b.word_len-2] + 2.0*b.radius*abs(1.0/z) ) {
        if (b.word_len < test_ball.word_len) {
          stack.push_back( act_on_right(0, b) );
          stack.push_back( act_on_right(1, b) );
        }
        min_dists.push_back(d);
        //std::cout << "New beginning min dist of " << d << " with " << Bitword(b.word, b.word_len) << "\n";
      }
    }
  }
  return min_dists.back();
}
    
/**************************************************************************
 * given a ball and bitwords, subdivide into at most two balls, 
 * together with all the bitwords which agree with them
 * ************************************************************************/
std::vector< std::pair<Ball,std::vector<Bitword> > > ifs::subdivide_ball_with_bitwords(const std::pair<Ball,std::vector<Bitword> >& b) {
  std::vector< std::pair<Ball,std::vector<Bitword> > > ans(2);
  ans[0].first = act_on_right(0, b.first);
  ans[0].second.resize(0);
  ans[1].first = act_on_right(1, b.first);
  ans[1].second.resize(0);
  for (int i=0; i<(int)b.second.size(); ++i) {
    int digit = b.second[i].reverse_get(b.first.word_len);
    ans[digit].second.push_back(b.second[i]);
  }
  if (ans[1].second.size() == 0) ans.erase(ans.begin()+1);
  if (ans[0].second.size() == 0) ans.erase(ans.begin());
  return ans;
}
/**************************************************************************
 * compute the minimum distance between two sets of balls
 * ************************************************************************/
double ifs::distance_between_bitword_sets(const std::vector<Bitword>& B1, 
                                          const std::vector<Bitword>& B2) {
  std::cout << "Finding the distance between the bitword sets:\n";
  for (int i=0; i<(int)B1.size(); ++i) {
    std::cout << i << ": " << B1[i] << "\n";
  }
  for (int i=0; i<(int)B2.size(); ++i) {
    std::cout << i << ": " << B2[i] << "\n";
  }
  double min_r;
  minimal_enclosing_radius(min_r);
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  Ball f = act_on_left(0, initial_ball);
  Ball g = act_on_left(1, initial_ball);
  std::vector<std::pair<Ball, std::vector<Bitword> > > B1_initials(0);
  B1_initials.push_back( std::make_pair( f, std::vector<Bitword>() ) );
  B1_initials.push_back( std::make_pair( g, std::vector<Bitword>() ) );
  std::vector<std::pair<Ball, std::vector<Bitword> > > B2_initials(0);
  B2_initials.push_back( std::make_pair( f, std::vector<Bitword>() ) );
  B2_initials.push_back( std::make_pair( g, std::vector<Bitword>() ) );
  for (int i=0; i<(int)B1.size(); ++i) {
    int digit = B1[i].reverse_get(0);
    B1_initials[digit].second.push_back( B1[i] );
  }
  if (B1_initials[1].second.size() == 0) B1_initials.erase( B1_initials.begin() + 1 );
  if (B1_initials[0].second.size() == 0) B1_initials.erase( B1_initials.begin() );
  for (int i=0; i<(int)B2.size(); ++i) {
    int digit = B2[i].reverse_get(0);
    B2_initials[digit].second.push_back( B2[i] );
  }
  if (B2_initials[1].second.size() == 0) B2_initials.erase( B2_initials.begin() + 1 );
  if (B2_initials[0].second.size() == 0) B2_initials.erase( B2_initials.begin() );
  //now the initials contain all the initial possibilities
  //build all the pairs
  std::deque< std::pair< std::pair<Ball,std::vector<Bitword> >, std::pair<Ball,std::vector<Bitword> > > > stack(0);
  for (int i=0; i<(int)B1_initials.size(); ++i) {
    for (int j=0; j<(int)B2_initials.size(); ++j) {
      stack.push_back( std::make_pair( B1_initials[i], B2_initials[j] ) );
    }
  }
  //now recursively search through
  std::vector<double> min_dists(0);
  while (stack.size() > 0) {
    std::pair< std::pair<Ball,std::vector<Bitword> >, std::pair<Ball,std::vector<Bitword> > > b = stack.back();
    stack.pop_back();
    
    std::cout << "Considering the ball pair: " << b.first.first << " " << b.second.first << " with " 
                                               << b.first.second.size() << " and " << b.second.second.size() << " bitwords\n";
    
    double d = abs(b.first.first.center - b.second.first.center) - 2.0*b.first.first.radius;
    int len = b.first.first.word_len;

    bool full_depth = (len == b.second.second[0].len);
    bool potentially_good_distance = ((int)min_dists.size() == 0) || 
                                     ((int)min_dists.size() < len && d < min_dists.back() + 4*b.first.first.radius*abs(1.0/z));
    if (potentially_good_distance) {
      //update the recorded distances
      if ((int)min_dists.size() < len) min_dists.push_back(d);
      else if (d < min_dists.back()) min_dists.back() = d;
      
      if (full_depth) continue;
      
      //actually do the subdivision 
      std::vector< std::pair<Ball,std::vector<Bitword> > > b1_subs = subdivide_ball_with_bitwords(b.first);
      std::vector< std::pair<Ball,std::vector<Bitword> > > b2_subs = subdivide_ball_with_bitwords(b.second);
      for (int i=0; i<(int)b1_subs.size(); ++i) {
        for (int j=0; j<(int)b2_subs.size(); ++j) {
          stack.push_front(std::make_pair( b1_subs[i], b2_subs[j] ) );
        }
      }
    }
  }
  return min_dists.back();
}
    
/**************************************************************************
 * check that all the balls have the given prefix
 **************************************************************************/
bool check_ball_prefixes(std::vector<Ball>& balls, 
                         std::vector<int>& L,
                         Bitword& b) {
  for (int m=0; m<(int)L.size(); ++m) {
    if (Bitword(balls[L[m]].word,
                balls[L[m]].word_len).prefix(b.len) != b) {
      return false;
    }
  }
  return true;
}


/*************************************************************************
 * at each depth up through n_depth, find the closest pair of balls 
 * (across fL,gL) which don't intersect
 *************************************************************************/
std::vector<double> ifs::shortest_nonintersection_distances(int n_depth, 
                                                            Ball& initial_ball, 
                                                            int verbose) {
  //at index i-1, we store the minimum distance at depth i
  std::vector<double> ans(n_depth, -1);
  Ball bf = act_on_right(0, initial_ball);
  Ball bg = act_on_right(1, initial_ball);
  std::deque<std::pair<Ball,Ball> > stack(0);
  stack.push_back(std::make_pair(bf,bg));
  while (stack.size() > 0) {
    std::pair<Ball,Ball> bp = stack.back();
    stack.pop_back();
    double d = (abs(bp.first.center - bp.second.center) - bp.first.radius) - bp.second.radius;
    int b_depth = bp.first.word_len;
    bool subdivide = (b_depth < n_depth) && (ans[b_depth-1] < 0 || d < ans[b_depth-1] + 4*bp.first.radius);
    bool new_min = (d > 0) && (ans[b_depth-1] < 0 || d < ans[b_depth-1]);
    if (new_min) {
      ans[b_depth-1] = d;
      if (verbose>0) {
        std::cout << "New min at depth " << b_depth << " of " << d << " with " 
                  << Bitword(bp.first.word, bp.first.word_len) << " " << Bitword(bp.second.word, bp.second.word_len) << "\n";
      }
    }
    if (subdivide) {
      Ball b1f = act_on_right(0, bp.first);
      Ball b1g = act_on_right(1, bp.first);
      Ball b2f = act_on_right(0, bp.second);
      Ball b2g = act_on_right(1, bp.second);
      stack.push_front(std::make_pair(b1f, b2f));
      stack.push_front(std::make_pair(b1f, b2g));
      stack.push_front(std::make_pair(b1g, b2f));
      stack.push_front(std::make_pair(b1g, b2g));
    }
  }
  return ans;
}










/**************************************************************************
 * functions to certify that the linear semiconjugacy is a conjugacy
 **************************************************************************/
bool ifs::certify_linear_conjugacy(double& epsilon, int n_depth, bool rigorous, int verbose) {

  if (abs(z) > 1.0/sqrt(2.0)) return false;
  if (z.real() < 0 || z.imag() < 0) return false;
  if (arg(z) < 0.15) return false;

  ifs temp_IFS;
  temp_IFS.set_params(z,z);
  temp_IFS.depth = n_depth;
  
  double min_r;
  if (!temp_IFS.minimal_enclosing_radius(min_r)) return false;
  if (!temp_IFS.circ_connected(min_r)) return false;
  
  if (verbose>0) {
    std::cout << "Certifying linear conjugacy at " << z << "\n";
  }
  
  //compute the balls, plus the pairs of balls which touch
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  std::vector<Ball> balls(0);
  temp_IFS.compute_balls(balls, initial_ball, n_depth);
  
  //compute the intersection pairs
  std::vector<std::pair<Ball,Ball> > intersection_pairs = temp_IFS.compute_intersection_pairs(n_depth, initial_ball, true);
  
  if (verbose>0) {
    std::cout << "Intersection balls:\n";
    for (int i=0; i<(int)intersection_pairs.size(); ++i) {
      std::cout << i << ": " << intersection_pairs[i].first << ",\n   " << intersection_pairs[i].second << "\n";
    }
  }
  
  //a 00 pair means a pair in which one of the balls has duplicate first letters
  //a 01 pair is anything else
  bool found_00_pair = false;
  bool found_01_pair = false;
  for (int i=0; i<(int)intersection_pairs.size(); ++i) {
    std::pair<Ball,Ball>& p = intersection_pairs[i];
    if (p.first.word_len < n_depth) continue;
    //these bools are true if the ball has duplicate first letters
    bool b1 = (p.first.word[p.first.word_len-1] == p.first.word[p.first.word_len-2]);
    bool b2 = (p.second.word[p.second.word_len-1] == p.second.word[p.second.word_len-2]);
    if (b1 || b2) {
      found_00_pair = true;
    } else {
      found_01_pair = true;
    }
  }
  
  if (verbose>0) std::cout << "found_00_pair: " << found_00_pair << "\nfound_01_pair: " << found_01_pair << "\n"; 
  
  //if there are both pairs with at least one ball with a 00/11 prefix
  //and pairs with 01/10, then we need to give up
  if (found_00_pair && found_01_pair) return false;
  
  //If the only intersecting pairs have at least one of the balls
  //with a 00 or 11 prefix, then we know we're conjugate since this is 
  //the case with one discontinuity.  The epsilon in this case 
  //is the one which prevents any 01/10 pairs from touching
  if (found_00_pair) {
    //get the closest 01 distance
    double d = temp_IFS.nonduplicate_first_letter_distance(n_depth, initial_ball, 0);
    //subtract the radius of the balls
    d -= 2*min_r*pow(abs(z), n_depth);
    
    //a bound on the derivative of the distance between centers is 
    //2 (1/Sqrt[2])/(1/Sqrt[2] - 1)^2 ~= 16.48528
    //the radii of the disks changes too; it can change as 
    //fast as 11.65685
    double deriv = 16.4853 + 2*11.6569;
    
    //compute epsilon from the distance
    epsilon = d/deriv;
    
    if (verbose>0) std::cout << "One discontinuity case; d = " << d << "; epsilon = " << epsilon << "\n";
    
    return true;
  }
  
  //if all the pairs are not 00/11, then we can go
  
  if (verbose>0) {
    std::cout << "Two discontinuities case\n";
  }
  
  //get the boundary -- normally, it is guaranteed to start at the g->f junction,
  //but we only need the section which runs from the rightmost 1-ball to 
  //the 0 balls, so we ask the function for that interval
  std::vector<int> boundary_indices(0);
  std::vector<Ball> boundary_balls(0);
  std::vector<Bitword> boundary_bitwords(0);
  non_grid_ball_boundary_indices(boundary_indices, balls, intersection_pairs, true, verbose-1);
  boundary_balls.resize(boundary_indices.size());
  boundary_bitwords.resize(boundary_indices.size());
  for (int i=0; i<(int)boundary_indices.size(); ++i) {
    boundary_balls[i] = balls[boundary_indices[i]];
    boundary_bitwords[i] = Bitword(boundary_balls[i].word, boundary_balls[i].word_len);
  }
  
  if (verbose>0) {
    std::cout << "Boundary balls:\n";
    for (int i=0; i<(int)boundary_balls.size(); ++i) {
      std::cout << i << " (" << boundary_indices[i] << "): " << boundary_balls[i] << "\n";
    }
  }
  
  //get the convex hull
  //std::vector<int> ch;
  //std::vector<cpx> ch_points;
  //std::vector<halfspace> ch_halfspaces;
  //ball_convex_hull(ch, ch_points, ch_halfspaces, boundary_balls);
  
  //if (verbose>0) {
  //  std::cout << "Convex hull:\n";
  //  for (int i=0; i<(int)ch.size(); ++i) {
  //    std::cout << i << ": " << ch[i] << "\n";
  //  }
  //}
  
  //since we have restricted the argument of the parameter, we only need the 
  //following: cyclically in order x,y,z,w,a, where x,y,z,w,a are certified 
  //prefixes of points on the boundary and:
  // x = 1^k*
  // y = W e^m (1-e)*
  // z = W e^(m+1) *
  // w = W e^m (1-e)*
  // a = 0*
  //where W begins with fewer than k 1's, and e = 0 or 1
  //the entire range along y,z,w must be certified to have the prefix We
  
  //first find the location of the largest number of 
  //1's -- i.e. a ball with prefix 1^k such that every ball it touches 
  //also has prefix 1^k; it'll probably be correct just to 
  //find the ball with the largest prefix of 1's and use that
  int max_num_1s_i = 0;
  int max_num_1s = 0;
  for (int i=1; i<(int)boundary_bitwords.size(); ++i) {
    int num_1s = boundary_bitwords[i].constant_prefix_size(1);
    if (num_1s > max_num_1s) {
      max_num_1s = num_1s;
      max_num_1s_i = i;
    }
  }
  //get all the balls which intersect the maximal-1 ball
  std::vector<int> max_num_1s_intersectors = balls_which_intersect_ball(boundary_balls[max_num_1s_i], 
                                                                        intersection_pairs, 
                                                                        verbose);
  //find the minimal size of the 1s prefix of these balls
  int k=max_num_1s;
  for (int i=0; i<(int)max_num_1s_intersectors.size(); ++i) {
    Bitword b(balls[max_num_1s_intersectors[i]].word, 
              balls[max_num_1s_intersectors[i]].word_len);
    int this_k = b.constant_prefix_size(1);
    if (this_k < k) {
      k = this_k;
    }
  }
  
  if (verbose>0) {
    std::cout << "Found maximal 1s prefix: " << boundary_balls[max_num_1s_i] << "; k = " << k << "\n";
  }
  
  //now we need to scan through, and for each prefix W, try to find three patterns 
  //y,z,w as desired 
  
  //find the beginning of the search interval; this will be the first 
  //ball which doesn't touch anything with the prefix 1^k
  int search_interval_start_i = max_num_1s_i;
  while (true) {
    std::vector<int> intersecting_balls = balls_which_intersect_ball(boundary_balls[search_interval_start_i],
                                                                     intersection_pairs, 
                                                                     verbose);
    bool ok=true;
    for (int i=0; i<(int)intersecting_balls.size(); ++i) {
      if (Bitword(balls[intersecting_balls[i]].word, 
                  balls[intersecting_balls[i]].word_len).constant_prefix_size(1) >= k) {
        ok = false;
        break;
      }
    }
    if (ok) break;
    if (search_interval_start_i == (int)boundary_bitwords.size()-1) return false;
    ++search_interval_start_i;
  }
  
  if (verbose>0) {
    std::cout << "Found search start index: " << search_interval_start_i << " (in balls: " << boundary_indices[search_interval_start_i] << ")\n";
  }
  
  int current_prefix_length = 2;
  while (true) {
    if (verbose>0) {
      std::cout << "Prefix length: " << current_prefix_length << "\n";
    }
    
    //scan through the list and find all prefixes of the given length
    //the prefix includes everything of the form We (including the e)
    int j=search_interval_start_i;
    std::vector<Bitword> prefixes(0);
    std::vector<std::pair<int,int> > prefix_ranges(0);  //this is [start,stop)
    while (true) {
      //check if we are done
      if (j >= (int)boundary_bitwords.size()) break;
      //get the prefix
      prefixes.push_back(boundary_bitwords[j].prefix(current_prefix_length));
      prefix_ranges.push_back(std::make_pair(j, -1));
      //go until the prefix isn't the prefix
      while (j < (int)boundary_bitwords.size() && boundary_bitwords[j].prefix(current_prefix_length) == prefixes.back()) ++j;
      prefix_ranges.back().second = j;
    }
    if (verbose>0) {
      std::cout << "Found prefixes:\n";
      for (int i=0; i<(int)prefixes.size(); ++i) {
        std::cout << i << ": " << prefixes[i] << " range: " << prefix_ranges[i].first << " " << prefix_ranges[i].second << "\n";
      }
    }
    //for each prefix, find the list of numbers of e's which occur
    //we need it to be <=m, then >m, then <=m again
    bool some_ok_prefix = false;
    for (int i=0; i<(int)prefixes.size(); ++i) {
      int prefix_range_length = prefix_ranges[i].second - prefix_ranges[i].first;
      //if the range length is at most 8, it cannot possibly work
      if (prefix_range_length <= 8) continue;
      some_ok_prefix = true;
      
      //get how many e's there are in each word; the bool indicates whether it's certified
      std::vector<std::pair<int,bool> > e_lengths(prefix_range_length);
      for (int j=0; j<(int)prefix_range_length; ++j) {
        int bd_ind = prefix_ranges[i].first+j;
        e_lengths[j].first = boundary_bitwords[bd_ind].continuance_of_prefix_last_letter(current_prefix_length);
        if (!rigorous) {
          e_lengths[j].second = true;
        } else {
          std::vector<int> intersecting_balls = balls_which_intersect_ball(boundary_balls[bd_ind],
                                                                         intersection_pairs,
                                                                         verbose);
          Bitword full_prefix = Bitword(boundary_balls[bd_ind].word,
                                        boundary_balls[bd_ind].word_len).prefix(current_prefix_length + e_lengths[j].first + 1);
          e_lengths[j].second = check_ball_prefixes(balls, intersecting_balls, full_prefix);          
        }
      }
      if (verbose>0) {
        std::cout << "For prefix " << i << ": " << prefixes[i] << "\n";
        std::cout << "e lengths:\n";
        for (int j=0; j<(int)e_lengths.size(); ++j) {
          std::cout << e_lengths[j].first << (e_lengths[j].second ? "* " : " ");
        }
        std::cout << "\n";
      }
      
      //find the peak
      int p1, p2, p3;
      find_integer_peak(e_lengths, p1, p2, p3);
      if (verbose>0) std::cout << "Found the peak: " << p1 << " " << p2 << " " << p3 << "\n";
      if (p1 < 0) {
        continue;
      } 
      if (!rigorous) {
        epsilon = 0;
        return true;
      }
      
      //we now have a rigorous conjugacy
      //we just need to determine the epsilon
      //it is sufficient to find a ball such that the intersection
      //pairs are smaller or remain the same
      std::vector<double>  ND = shortest_nonintersection_distances(n_depth, initial_ball, 0);
      
      if (verbose>0) {
        std::cout << "Nonintersection distances:\n";
        for (int j=0; j<(int)ND.size(); ++j) {
          std::cout << j+1 << ": " << ND[j] << "\n";
        }
      }
      
      //find the bound on the distance for z: for each ball distance, 
      //we get derivative estimate, and we know that the distance 
      //can be covered by that derivative plus the rate at which the radii 
      //can change (times 2), which is 11.6569
      double min_z_epsilon = -1;
      double Z = 1.0/sqrt(2.0);
      double Zm1 = Z-1.0;
      double radius_deriv = 11.6569;
      for (int j=0; j<(int)ND.size(); ++j) {
        double diff_deriv = (2.0 + 2.0*(-1.0 + (j+1)*Zm1)*pow(Z,j+1)) / (Zm1*Zm1);
        double total_deriv = diff_deriv + 2*radius_deriv;
        double z_epsilon = ND[j] / total_deriv;
        if (verbose>0) std::cout << "Computed total deriv " << total_deriv << " and epsilon from depth " << j+1 << ": " << z_epsilon << "\n";
        if (z_epsilon > 0 && (min_z_epsilon < 0 || z_epsilon < min_z_epsilon)) {
          min_z_epsilon = z_epsilon;
        }
      }
      
      epsilon = min_z_epsilon;      
      
      return true;
      
      
    }
    if (!some_ok_prefix) break;
    ++current_prefix_length;
  }

  return false;
}





/****************************************************************************
 * find coordinates theta,lambda using the kneading invariant (address) 
 * of the points of intersection
 * *************************************************************************/
bool ifs::coordinates_from_kneading(double& theta, double& lambda, 
                                    int n_depth, int verbose) {
  //do sanity checks
  if (abs(z) > 1.0/sqrt(2.0)) return false;
  if (z.real() < 0 || z.imag() < 0) return false;
  if (arg(z) < 0.15) return false;

  ifs temp_IFS;
  temp_IFS.set_params(z,z);
  temp_IFS.depth = n_depth;
  
  double min_r;
  if (!temp_IFS.minimal_enclosing_radius(min_r)) return false;
  if (!temp_IFS.circ_connected(min_r)) return false;
  
  if (verbose>0) {
    std::cout << "Computing coordinates at " << z << "\n";
  }
  
  //compute the balls, plus the pairs of balls which touch
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  std::vector<Ball> balls(0);
  temp_IFS.compute_balls(balls, initial_ball, n_depth);
  
  //compute the intersection pairs
  std::vector<std::pair<Ball,Ball> > intersection_pairs = temp_IFS.compute_intersection_pairs(n_depth, initial_ball, true);
  
  if (verbose>1) {
    std::cout << "Intersection balls:\n";
    for (int i=0; i<(int)intersection_pairs.size(); ++i) {
      std::cout << i << ": " << intersection_pairs[i].first << ",\n   " << intersection_pairs[i].second << "\n";
    }
  }
  
  //the only thing we care about the the address of the balls
  //as we go from g to f, so we get that pair
  Bitword g_bitword, f_bitword;
  non_grid_ball_boundary_transition_indices(g_bitword, f_bitword,
                                            balls, intersection_pairs, verbose-1);
  if (verbose>0) {
    std::cout << "Got transition pair: " << g_bitword << " " << f_bitword << "\n";
  }

  theta = -1;
  lambda = -1;
  
  return true;
}


/****************************************************************************
 * find the dynamical lamination
 * optionally, find only the addresses of (leaves hiding) the endpoints
 * of the leaf ell_f
 * **************************************************************************/
bool ifs::dynamical_lamination(int n_depth, 
                               bool only_hiding_ell_f, 
                               int verbose) {
  return true;
}








