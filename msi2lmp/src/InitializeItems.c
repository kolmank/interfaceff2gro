/*
*   This function fills in the keyword field, the number of members for each
*   item and the number of parameters for each item
*
*/

#include "msi2lmp.h"
#include "Forcefield.h"

#include <string.h>

void InitializeItems(void)
{
  /* ATOM TYPES */
  strcpy_s(ff_atomtypes.keyword,sizeof(ff_atomtypes.keyword),"#atom_types");
  ff_atomtypes.number_of_members = 1;
  ff_atomtypes.number_of_parameters = 1;

  /* EQUIVALENCE */

  strcpy_s(equivalence.keyword,sizeof(equivalence.keyword),"#equivalence");
  equivalence.number_of_members = 6;
  equivalence.number_of_parameters = 0;

  /* NON-BOND */

  strcpy_s(ff_vdw.keyword,sizeof(ff_vdw.keyword),"#nonbond");
  ff_vdw.number_of_members = 1;
  ff_vdw.number_of_parameters = 2;

  /* BOND */

  ff_bond.number_of_members = 2;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy_s(ff_bond.keyword,sizeof(ff_bond.keyword),"#quadratic_bond");
    ff_bond.number_of_parameters = 2;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy_s(ff_bond.keyword,sizeof(ff_bond.keyword),"#quartic_bond");
    ff_bond.number_of_parameters = 4;
  }

  /* MORSE */

  if (forcefield & FF_TYPE_CLASS1) {
    ff_morse.number_of_members = 2;
    strcpy_s(ff_morse.keyword,sizeof(ff_morse.keyword),"#morse_bond");
    ff_morse.number_of_parameters = 3;
  }

  /* ANGLE */

  ff_ang.number_of_members = 3;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy_s(ff_ang.keyword,sizeof(ff_ang.keyword),"#quadratic_angle");
    ff_ang.number_of_parameters = 2;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy_s(ff_ang.keyword,sizeof(ff_ang.keyword),"#quartic_angle");
    ff_ang.number_of_parameters = 4;
  }

  /* TORSION */

  ff_tor.number_of_members = 4;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy_s(ff_tor.keyword,sizeof(ff_tor.keyword),"#torsion_1");
    ff_tor.number_of_parameters = 3;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy_s(ff_tor.keyword,sizeof(ff_tor.keyword),"#torsion_3");
    ff_tor.number_of_parameters = 6;
  }

  /* OOP */

  ff_oop.number_of_members = 4;
  if (forcefield & (FF_TYPE_CLASS1|FF_TYPE_OPLSAA)) {
    strcpy_s(ff_oop.keyword,sizeof(ff_oop.keyword),"#out_of_plane");
    ff_oop.number_of_parameters = 3;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    strcpy_s(ff_oop.keyword,sizeof(ff_oop.keyword),"#wilson_out_of_plane");
    ff_oop.number_of_parameters = 2;
  }

  if (forcefield & FF_TYPE_CLASS2) {
    /* BOND-BOND */

    strcpy_s(ff_bonbon.keyword,sizeof(ff_bonbon.keyword),"#bond-bond");
    ff_bonbon.number_of_members = 3;
    ff_bonbon.number_of_parameters = 1;

    /* BOND-ANGLE */

    strcpy_s(ff_bonang.keyword,sizeof(ff_bonang.keyword),"#bond-angle");
    ff_bonang.number_of_members = 3;
    ff_bonang.number_of_parameters = 2;

    /* ANGLE-TORSION */

    strcpy_s(ff_angtor.keyword,sizeof(ff_angtor.keyword),"#angle-torsion_3");
    ff_angtor.number_of_members = 4;
    ff_angtor.number_of_parameters = 6;

    /* ANGLE-ANGLE-TORSION */

    strcpy_s(ff_angangtor.keyword,sizeof(ff_angangtor.keyword),"#angle-angle-torsion_1");
    ff_angangtor.number_of_members = 4;
    ff_angangtor.number_of_parameters = 1;

    /* END-BOND-TORSION */

    strcpy_s(ff_endbontor.keyword,sizeof(ff_endbontor.keyword),"#end_bond-torsion_3");
    ff_endbontor.number_of_members = 4;
    ff_endbontor.number_of_parameters = 6;

    /* MID-BOND-TORSION */

    strcpy_s(ff_midbontor.keyword,sizeof(ff_midbontor.keyword),"#middle_bond-torsion_3");
    ff_midbontor.number_of_members = 4;
    ff_midbontor.number_of_parameters = 3;

    /* ANGLE-ANGLE */

    strcpy_s(ff_angang.keyword,sizeof(ff_angang.keyword),"#angle-angle");
    ff_angang.number_of_members = 4;
    ff_angang.number_of_parameters = 1;

    /* BOND-BOND-1-3 */

    strcpy_s(ff_bonbon13.keyword,sizeof(ff_bonbon13.keyword),"#bond-bond_1_3");
    ff_bonbon13.number_of_members = 4;
    ff_bonbon13.number_of_parameters = 1;
  }
}
