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
  return 0;
}







/****************************************************************************
* compute boundary balls not using a trap grid
****************************************************************************/
void non_grid_ball_boundary_indices(std::vector<int>& boundary, 
                                    std::vector<Ball>& balls, 
                                    std::vector<std::pair<Ball,Ball> >& intersection_pairs,
                                    int verbose) {
  boundary.resize(0);
  
  //find the ball with the greatest x coord 
  int max_x_ball_i = 0;
  for (int i=0; i<(int)balls.size(); ++i) {
    if (balls[i].center.real() > balls[max_x_ball_i].center.real()) {
      max_x_ball_i = balls[i].center.real();
    }
  }
  
  //so balls[max_x_ball_i] must have angle 0 exposed
  int current_ball_i = max_x_ball_i;
  std::vector<int> intersecting_balls = balls_which_intersect_ball(balls[current_ball_i], intersection_pairs);
  //get the ball which is previous in the boundary
  int prev_ball_i = last_intersecting_ball_before_angle(balls, current_ball_i, intersecting_balls, 0.0);  
  int next_ball_i;
  
  while (true) {
    //get the next ball after the previous one
    next_ball_i = next_intersecting_ball_after_ball(balls, current_ball_i, prev_ball_i, intersecting_balls);
    
    //we might have looped around
    if (boundary.size() > 2 &&
        current_ball_i == boundary[0] &&
        next_ball_i == boundary[1]) {
      break;
    }
    
    //if we haven't looped, we should put on the current ball
    boundary.push_back(current_ball_i);
    
    //next becomes current and current becomes previous
    prev_ball_i = current_ball_i;
    current_ball_i = next_ball_i;
    
    //compute the intersecting balls
    intersecting_balls = balls_which_intersect_ball(balls[current_ball_i], intersection_pairs);
  }
  
}


/**************************************************************************
 * functions to certify that the linear semiconjugacy is a conjugacy
 **************************************************************************/
bool ifs::certify_linear_conjugacy(double& epsilon, int n_depth, int verbose) {

  if (abs(z) > 1.0/sqrt(2.0)) return false;
  if (z.real() < 0 || z.imag() < 0) return false;

  ifs temp_IFS;
  temp_IFS.set_params(z,z);
  temp_IFS.depth = n_depth;
  
  double min_r;
  if (!temp_IFS.minimal_enclosing_radius(min_r)) return false;
  if (!temp_IFS.circ_connected(min_r)) return false;
  
  //compute the balls, plus the pairs of balls which touch
  Ball initial_ball(0.5,(z-1.0)/2.0,(1.0-w)/2.0,min_r);
  std::vector<Ball> balls(0);
  temp_IFS.compute_balls(balls, initial_ball, n_depth);
  
  //compute the intersection pairs
  std::vector<std::pair<Ball,Ball> > intersection_pairs = temp_IFS.compute_intersection_pairs(n_depth, initial_ball, true);
  
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
    double deriv = 16.4853 + 11.6569;
    
    //compute epsilon from the distance
    epsilon = d/deriv;
    
    return true;
  }
  
  //if all the pairs are not 00/11, then we can go
  
  //get the boundary (it is guaranteed to start at the g->f junction)
  std::vector<int> boundary_indices(0);
  std::vector<Ball> boundary_balls(0);
  non_grid_ball_boundary_indices(boundary_indices, balls, intersection_pairs);
  boundary_balls.resize(boundary_indices.size());
  for (int i=0; i<(int)boundary_indices.size(); ++i) {
    boundary_balls[i] = balls[boundary_indices[i]];
  }
  
  //get the convex hull
  std::vector<int> ch;
  std::vector<cpx> ch_points;
  std::vector<halfspace> ch_halfspaces;
  ball_convex_hull(ch, ch_points, ch_halfspaces, boundary_balls);
  
  
  
  
  

  return true;
}























