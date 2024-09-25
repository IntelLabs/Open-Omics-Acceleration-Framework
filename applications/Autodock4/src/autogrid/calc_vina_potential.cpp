/*
$ID$
Copyright (C) 2011 The Scripps Research Institute. All rights reserved.

 AutoGrid is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <math.h>
#include <string.h>
#include "autogrid.h"
#include "atom_parameter_manager.h"
#include "calc_vina_potential.h"


int calc_vina_potential( const int num_ligand_types, 
                          const char ligand_types[][3],
                          const int num_receptor_types, 
                          const char receptor_types[][3],
          /* not const */ EnergyTables *const p_et) 
{
    // loop up to a maximum distance of  (NEINT * INV_A_DIV), 
    //                          usually    2048 * 0.01,       or 20.48 Angstroms
    //int debug=1; 
    //from autodock_vina_1_1_1/src/main.cpp,line 394
    float wt_gauss1 = -0.035579; //vina:everything.cpp line321
    float wt_gauss2 = -0.005156; //vina:everything.cpp line343
    float wt_repulsion = 0.840245;
    float wt_hydrogen = -0.587439;
    float wt_hydrophobic = -0.035069; //C_H,F_H,Cl_H,Br_H,I_H
    // xs_vdw_radii from atom_constants.h now in read_parameter_library.cc
    //float C_H = 1.9;//C_P
    //float N_P = 1.8;//N_D,N_A,N_DA
    //float O_P = 1.7;//O_D,O_A,O_DA
    //float S_P = 2.0;
    //float P_P = 2.1;
    //float F_H = 1.5;
    //float Cl_H = 1.8;
    //float Br_H = 2.0;
    //float I_H = 2.2;
    //float Met_D = 1.2; //metal_donor:Mg,Mn,Zn,Ca,Fe,Cl,Br
    //float Met_non_ad = 1.75;//metal_non_ad:Cu,Fe,Na,K,Hg,Co,U,Cd,Ni
    double rddist; 
    double delta_e = 0.0;
    //variables added Aug10 2011
    int ia, i,   indx_r;
    double r, Rij, map_Rij;

    //                                                               ia_dist
    //    interatom_distance                                   |.................|
    //    interatom_distance - xs_radius(t1) -xs_radius(t2)    .--|........|-----.  
    //                                                       at1               at2
    //vina distance from current gridpt to atom ia:        xs_rad1  rddist   xs_rad2
    //0. process receptor to setup type[ia], coords[ia], xs_rad[ia]
    //1. setup map types
    //2. setup the energy_lookup tables
    //3. loop over all the maps 
    //4.     loop over all pts in current map_ia
    //5.        loop over all the receptor atoms adding to this pt
    //             rdist: interatom_distance from atom coords to current grid pt
    //                rddist based on types: ia_dist - (xs_rad1 + xs_rad2)
    // terms:
    //     e_attractive:
    //         delta_e = rgauss1*exp(-((rddist)/0.5)**2) + rgauss2*exp(-((rddist-3.)/2.)**2);
    //         energy_lookup[i][indx_r][ia] += delta_e
    //     e_repulsive:
    //         if (rddist<0.0){
    //             delta_e = rrepulsive*(rddist**2);
    //             energy_lookup[i][indx_r][ia] += delta_e;
    //         };
    //     e_hbond:
    //         (1)set ihb from types; it is set to 1 if pair of types is suitable for hbond int.
    //         ihb = 0;
    //         if (ihb>0){
    //             if (rddist<0.7)
    //                  delta_e = 1*weight_hydrogen;
    //                  energy_lookup[i][indx_r][ia] += delta_e;
    //             if ((-0.7<rddist) && (rddist<0.))
    //                  delta_e =(rddist/0.7)*weight_hydrogen;
    //                  energy_lookup[i][indx_r][ia] -= delta_e;
    //         }
    //     e_hydrophobic:
    //         //energy_lookup[atom_type[ia]][indx_r][map_ia]+= e_hphob
    //         (1)TODO: set ihb from types; it is set to 1 if pair of types is suitable for hydrophobic int.
    //         ihb = 0
    //         if (rddist<0.5){
    //              energy_lookup[i][indx_r][ia]+= 1*weight_hydrophobic;}
    //         else if ((0.5<rddist)&& (rdist<1.5)) {
    //              energy_lookup[i][indx_r][ia]+=(0.5-rddist)*weight_hydrophobic;
    //         };
    // to use in filling out the gridmap
    //????      gridmap[map_index].energy += energy_lookup[i][indx_r][i];
    //6.     output this map
    //
    /*canned receptor atom type number*/
    int hydrogen=-1; 
    int carbon=-1; 
    int arom_carbon=-1; 
    int oxygen=-1;
    int nitrogen=-1; 
    int fluorine=-1;
    int chlorine=-1;
    int bromine=-1; 
    int iodine=-1;
    int lig_hbond=-1;
    int rec_hbond=-1;
    ParameterEntry * lig_parm; 
    ParameterEntry * rec_parm; 

    for (ia=0; ia<num_ligand_types; ia++){
        lig_parm = apm_find(ligand_types[ia]);
        if ( lig_parm != NULL ) {
            //(void) fprintf ( logFile, "DEBUG: lig_parm->lig_index = %d", ia );
            (void) printf ( "DEBUG: lig_parm->lig_index = %d\n", ia );
        }
        for (i = 0;  i < num_receptor_types;  i++) {
            /*for each receptor_type*/
            //xA = gridmap[ia].xA[i];
            //xB = gridmap[ia].xB[i];
            //Rij = gridmap[ia].nbp_r[i];
            //epsij = gridmap[ia].nbp_eps[i];
            //ParameterEntry * rec_parm = apm_find(receptor_types[i]);
            rec_parm = apm_find(receptor_types[i]); //@@
            if (strcmp(receptor_types[i],"H")==0) { hydrogen=i; };
            if (strcmp(receptor_types[i],"C")==0) { carbon=i; };
            if (strcmp(receptor_types[i],"A")==0) { arom_carbon=i; };
            if (strcmp(receptor_types[i],"O")==0) { oxygen=i; };
            if (strcmp(receptor_types[i],"N")==0) { nitrogen=i; };
            if (strcmp(receptor_types[i],"F")==0) { fluorine=i; };
            if (strcmp(receptor_types[i],"Cl")==0) { chlorine=i; };
            if (strcmp(receptor_types[i],"Br")==0) { bromine=i; };
            if (strcmp(receptor_types[i],"I")==0) { iodine=i; };
            // get xs_radius for this probe type from parameter_library
            Rij = lig_parm->xs_radius; //see read_parameter_library
            lig_hbond = lig_parm->hbond;
            map_Rij = rec_parm->xs_radius; //see read_parameter_library
            rec_hbond = rec_parm->hbond;
            //@@TODO@@: add SER-OG,THR-OG, TYR_OH: X(1.2) Cl_H(1.8),Br_H(2.0),I_H(2.2),Met_D(1.2)
            /* loop over distance index, indx_r, from 0 to NDIEL */ /* GPF_MAP */
#ifdef DEBUG
            printf("%d-%d-building  Rij=%6.3lf, map_Rij=%10.8f for %s %s\n",ia,i, Rij, map_Rij, receptor_types[i], ligand_types[ia]);
#endif
            //(void) fprintf( logFile, "Calculating vina energies for %s-%s interactions (%d, %d).\n", gridmap[ia].type, receptor_types[i], ia, i );
            for (indx_r = 1;  indx_r < NDIEL;  indx_r++) {
                r  = angstrom(indx_r);
                // compute rddist:                                   map_Rij  rddist   Rij
                //  interatom_distance - xs_radius(t1) -xs_radius(t2)    .--|........|-----.  
                rddist =  r - (map_Rij + Rij);
                //use rddist for computing the vina component energies
                //@@TODO@@: replace with functions from vina..
                // NB: everything.cpp line 321: offset=0.0;width=0.5;cutoff=8
                //attraction:
                //8/2011: attraction from vina: everything.cpp lines 26-28
                //inline fl gaussian(fl x, fl width) {  
                //    return std::exp(-sqr(x/width));
                //}
                //8/2011: optimal_distance from vina: everything.cpp lines 94-96
                //inline fl optimal_distance(sz xs_t1, sz xs_t2) {
                //    return xs_radius(xs_t1) + xs_radius(xs_t2);
                //}
                //8/2011: vina: everything.cpp lines 98~107
                //struct gauss : public usable {
                //    fl offset; // added to optimal distance
                //    fl width;
                //    gauss(fl offset_, fl width_, fl cutoff_) : usable(cutoff_), offset(offset_), width(width_) {
                //        name = std::string("gauss(o=") + to_string(offset) + ", w=" + to_string(width) + ", c=" + to_string(cutoff) + ")";
                //    }
                //    fl eval(sz t1, sz t2, fl r) const {
                //        return gaussian(r - (optimal_distance(t1, t2) + offset), width);
                //    }
                //};
                // first gauss:  wt=-0.035579, offset=0., width=0.5 
                // second gauss: wt=-0.005156, offset=3., width=2.0
                delta_e = wt_gauss1 * exp(-pow(((rddist)/0.5),2)) + wt_gauss2 * exp(-pow(((rddist-3.)/2.),2));
                //at distance 'indx_r': interaction of receptor atomtype 'ia' - ligand atomtype 'i'
                //@@// 1: energy_lookup[i][indx_r][ia] += delta_e; 
                //8/2011: repulsion from vina: everything.cpp l 109-120
                //struct repulsion : public usable {
                //	fl offset; // added to vdw
                //	repulsion(fl offset_, fl cutoff_) : usable(cutoff_), offset(offset_) {
                //		name = std::string("repulsion(o=") + to_string(offset) + ")";
                //	}
                //	fl eval(sz t1, sz t2, fl r) const {
                //		fl d = r - (optimal_distance(t1, t2) + offset);
                //		if(d > 0) 
                //			return 0;
                //		return d*d;
                //	}
                //}                
                if (rddist<0){
                    delta_e = wt_repulsion*pow(rddist,2);
                    p_et->e_vdW_Hb[i][indx_r][ia] += delta_e;
                    //energy_lookup[i][indx_r][ia] += delta_e;
                }
                //hbond
                //add(1, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff // WEIGHT:  -0.587439
                // good_ = -0.7
                // bad_  = 0.0
                // cutoff  = 8
                // wt_hydrogen = -0.587439
                //struct non_dir_h_bond : public usable {
                //    fl good;
                //    fl bad;
                //    non_dir_h_bond(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
                //        name = std::string("non_dir_h_bond(g=") + to_string(good) + ", b=" + to_string(bad) + ")";
                //    }
                //    fl eval(sz t1, sz t2, fl r) const {
                //        if(xs_h_bond_possible(t1, t2))
                //            return slope_step(bad, good, r - optimal_distance(t1, t2));
                //        return 0;
                //    }
                //};


                if ((rec_hbond>2&& (lig_hbond==1||lig_hbond==2))||((rec_hbond==1||rec_hbond==2)&&lig_hbond>2)){ //check that types ia-i hbond
#ifdef DEBUG
                    printf(" processing gridmap= %d-hbonder i= %d\n",ia, i);
#endif
                    if (rddist<=0.7) { //what about EXACTLY 0.7?
                       delta_e = 1*wt_hydrogen;
                       p_et->e_vdW_Hb[i][indx_r][ia] += delta_e;
                       //energy_lookup[i][indx_r][ia] += delta_e;
                    }
                    if ((-0.7<rddist) && (rddist<=0.)){
                       delta_e =(rddist/0.7)*wt_hydrogen;
                       p_et->e_vdW_Hb[i][indx_r][ia] += delta_e;
                       //energy_lookup[i][indx_r][ia] -= delta_e;
                    }
                }
                // hydrophobic: check using index 'i' compared with 'carbon',
                //8/2011: hydrophobic from vina: everything.cpp l 134-145
                //struct hydrophobic : public usable {
                //fl good;
                //fl bad;
                //hydrophobic(fl good_, fl bad_, fl cutoff_) : usable(cutoff_), good(good_), bad(bad_) {
                //    name = "hydrophobic(g=" + to_string(good) + ", b=" + to_string(bad) + ", c=" + to_string(cutoff) + ")";
                //}
                //fl eval(sz t1, sz t2, fl r) const {
                //    if(xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
                //        return slope_step(bad, good, r - optimal_distance(t1, t2));
                //    else return 0;
                //}
                //};

                //	add(1, new hydrophobic(0.5, 1.5, cutoff)); // good, bad, cutoff // WEIGHT:  -0.035069
                // carbon/aromatic_carbon       to non-hbonder
                //@@TODO: add support for these other hydrophobic interactions:
                //if (((i==carbon)||(i==arom_carbon)||(i==fluorine)||(i==chlorine)||(i==bromine)||(i==iodine))
                //&& ((ia==carbon)||(ia==arom_carbon)||(ia==fluorine)||(ia==chlorine)||(ia==bromine)||(ia==iodine))) 
                //wt_hydrophobic=-0.035069
                if (((i==carbon)||(i==arom_carbon)) && ((ia==carbon)||(ia==arom_carbon)))
                 {
                    delta_e = 0.;
                    if (rddist<0.5) {
                       delta_e = 1*wt_hydrophobic;
                    } else if (rddist<1.5){
                       delta_e = (0.5-rddist)*wt_hydrophobic;
                    }
                    p_et->e_vdW_Hb[i][indx_r][ia] += delta_e;
                    //energy_lookup[i][indx_r][ia] += delta_e;
                } //i + ia both some type of carbon 
            } /*for each distance indx_r up to MAX_DIST*/ 
          } /*for each receptor_type*/ 
#ifdef DEBUG
            printf("END USE_VINA_POTENTIAL\n");
#endif
     } /* for each ligand type */
     printf("returning from USE_VINA_POTENTIAL\n");
     printf("last value=%f\n", p_et->e_vdW_Hb[i][indx_r][ia]);
     printf("last delta_e value=%f, rddist=%f, i=%d, indx_r=%d, ia=%d\n", delta_e, rddist,i,indx_r, ia);
     return 1;

//add(1, new num_tors_div()); // WEIGHT: 1.923 -- FIXME too close to limit?
} /* END calc_vinapotential*/
