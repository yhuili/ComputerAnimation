/*

  USC/Viterbi/Computer Science
  "Particle System" Assignment 3 starter code

*/


#ifndef _INPUT_H_
#define _INPUT_H_

/* Write a screenshot to the specified filename, in PPM format */
void saveScreenshot (int windowWidth, int windowHeight, char *filename);

// mouse & keyboard control
void mouseMotionDrag(int x, int y);
void mouseMotion (int x, int y);
void mouseButton(int button, int state, int x, int y);
void keyboardFunc (unsigned char key, int x, int y);
void SpecialKeys(int key, int x, int y);

#endif

