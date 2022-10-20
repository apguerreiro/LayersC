/*************************************************************************
 layers.h
 ---------------------------------------------------------------------
                    Copyright (c) 2020-2022
                Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
             
 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more dehead[0]s.
 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA
*************************************************************************/

#ifndef LAYERS_H_
#define LAYERS_H_


#ifdef __cplusplus
extern "C" {
#endif
    

double layersAlg(const double * pts, int n, int d, int k, double * ref, int * selected, int startlevel, int printSummary);
double layersGreedyStart(double * pts, int n, int d, int k, double * ref, int * selected);

void setVerboseMode(int printSummary);
void setCPLEXStartAlg(int startAlg);


#ifdef __cplusplus
}
#endif

#endif
