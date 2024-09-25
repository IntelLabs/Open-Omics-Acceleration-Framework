/*

 $Id: banner.cc,v 1.25 2014/06/12 01:44:07 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "banner.h"

/* local convenience function to return "OpenMP" or "" depending on C option */
static char* str_openmp(void); 

void banner( const char *const version_num, const int outlev, FILE *logFile )
{

/*----------------------------------------------------------------------------*/
/* Output banner...                                                           */
/*----------------------------------------------------------------------------*/

if(outlev>=LOGRUNV) {

    (void) fprintf(logFile, "      ________________________________________________________________\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "__________//___________________________/////___________________/____________\n");
    (void) fprintf(logFile, "_________/__/__________________________/____/__________________/____________\n");
    (void) fprintf(logFile, "________/____/___________/_____________/_____/_________________/____________\n");
    (void) fprintf(logFile, "________/____/__/_____/_/////___/////__/_____/__/////___/////__/___/________\n");
    (void) fprintf(logFile, "_______/______/_/_____/__/_____/_____/_/_____/_/_____/_/_____/_/_//_________\n");
    (void) fprintf(logFile, "_______////////_/_____/__/_____/_____/_/_____/_/_____/_/_______//_/_________\n");
    (void) fprintf(logFile, "_______/______/_/____//__/___/_/_____/_/____/__/_____/_/_____/_/___/________\n");
    (void) fprintf(logFile, "_______/______/__////_/___///___/////__/////____/////___/////__/____/_______\n\n");
    (void) fprintf(logFile, "      ________________________________________________________________\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "                                ______\n");
    (void) fprintf(logFile, "                               /      \\\n");
    (void) fprintf(logFile, "                              /        \\\n");
    (void) fprintf(logFile, "                             /          \\\n");
    (void) fprintf(logFile, "                             \\    /\\    /\n");
    (void) fprintf(logFile, "                              \\  /  \\  /\n");
    (void) fprintf(logFile, "                               \\/ /\\ \\/\n");
    (void) fprintf(logFile, "                                 /  \\\n");
    (void) fprintf(logFile, "                                /____\\\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "                  ______________________________________ \n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 |    AutoDock %-3.3s Release %-8s%s  |\n", version_num, version_num, str_openmp());
    (void) fprintf(logFile, "                 |            (C) 1989-2012             |\n");
    (void) fprintf(logFile, "                 |    The Scripps Research Institute    |\n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 |        Garrett M. Morris, TSRI       |\n");
    (void) fprintf(logFile, "                 |            Ruth Huey, TSRI           |\n");
    (void) fprintf(logFile, "                 |          Michael Pique, TSRI         |\n");
    (void) fprintf(logFile, "                 |        William E. Hart, Sandia       |\n");
    (void) fprintf(logFile, "                 |   R. Scott Halliday, Hewlett Packard |\n");
    (void) fprintf(logFile, "                 |        William Lindstrom, TSRI       |\n");
    (void) fprintf(logFile, "                 |           Max Chang, TSRI            |\n");
    (void) fprintf(logFile, "                 |        Alexander Gillet, TSRI        |\n");
    (void) fprintf(logFile, "                 |         Stefano Forli, TSRI          |\n");
    (void) fprintf(logFile, "                 |         Chenglong Li, OSU            |\n");
    (void) fprintf(logFile, "                 |         Huameng Li, OSU              |\n");
    (void) fprintf(logFile, "                 |        Richard K. Belew, UCSD        |\n");
    (void) fprintf(logFile, "                 |       David S. Goodsell, TSRI        |\n");
    (void) fprintf(logFile, "                 |        Arthur J. Olson, TSRI         |\n");
    (void) fprintf(logFile, "                 |______________________________________|\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "                  ______________________________________ \n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 | Automated Docking of Flexible Ligand |\n");
    (void) fprintf(logFile, "                 |  to Flexible Macromolecular Receptor |\n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 | For help, email %-19s |\n", PACKAGE_BUGREPORT);
    (void) fprintf(logFile, "                 |______________________________________|\n");
    (void) fprintf(logFile, "\n\n");
}
else {
    (void) fprintf(logFile, "          AutoDock %-3.3s Release %-8s%s\n", version_num, version_num, str_openmp());
    (void) fprintf(logFile, "         (C) 1989-2012 The Scripps Research Institute\n");
}

    (void) fprintf(logFile, "        AutoDock comes with ABSOLUTELY NO WARRANTY.\n");
// GNU BEGIN   (see maintenance script update_license_de-GNU)
    (void) fprintf(logFile, "        AutoDock is free software, and you are welcome\n");
    (void) fprintf(logFile, "        to redistribute it under certain conditions;\n");
    (void) fprintf(logFile, "        for details type 'autodock4 -C'\n\n");
// GNU END   (see maintenance script update_license_de-GNU)

}

static char* str_openmp(void) 
{ 
#ifdef _OPENMP
static char s[] = " [OpenMP]";
#else
static char s[] = "";
#endif
return s;
}

/*----------------------------------------------------------------------------*/

void show_copyright( FILE *const fp )
{

// GNU BEGIN   (see maintenance script update_license_de-GNU)
    (void) fprintf( fp, "GNU GENERAL PUBLIC LICENSE\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Version 2, June 1991\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Copyright (C) 1989, 1991 Free Software Foundation, Inc. 51 Franklin\n");
    (void) fprintf( fp, "Street, Fifth Floor, Boston, MA  02110-1301, USA\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Everyone is permitted to copy and distribute verbatim copies of this\n");
    (void) fprintf( fp, "license document, but changing it is not allowed.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "                       Preamble\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "The licenses for most software are designed to take away your freedom to\n");
    (void) fprintf( fp, "share and change it. By contrast, the GNU General Public License is\n");
    (void) fprintf( fp, "intended to guarantee your freedom to share and change free software--to\n");
    (void) fprintf( fp, "make sure the software is free for all its users. This General Public\n");
    (void) fprintf( fp, "License applies to most of the Free Software Foundation's software and\n");
    (void) fprintf( fp, "to any other program whose authors commit to using it. (Some other Free\n");
    (void) fprintf( fp, "Software Foundation software is covered by the GNU Lesser General Public\n");
    (void) fprintf( fp, "License instead.) You can apply it to your programs, too.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "When we speak of free software, we are referring to freedom, not price.\n");
    (void) fprintf( fp, "Our General Public Licenses are designed to make sure that you have the\n");
    (void) fprintf( fp, "freedom to distribute copies of free software (and charge for this\n");
    (void) fprintf( fp, "service if you wish), that you receive source code or can get it if you\n");
    (void) fprintf( fp, "want it, that you can change the software or use pieces of it in new\n");
    (void) fprintf( fp, "free programs; and that you know you can do these things.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "To protect your rights, we need to make restrictions that forbid anyone\n");
    (void) fprintf( fp, "to deny you these rights or to ask you to surrender the rights. These\n");
    (void) fprintf( fp, "restrictions translate to certain responsibilities for you if you\n");
    (void) fprintf( fp, "distribute copies of the software, or if you modify it.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "For example, if you distribute copies of such a program, whether gratis\n");
    (void) fprintf( fp, "or for a fee, you must give the recipients all the rights that you have.\n");
    (void) fprintf( fp, "You must make sure that they, too, receive or can get the source code.\n");
    (void) fprintf( fp, "And you must show them these terms so they know their rights.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "We protect your rights with two steps: (1) copyright the software, and\n");
    (void) fprintf( fp, "(2) offer you this license which gives you legal permission to copy,\n");
    (void) fprintf( fp, "distribute and/or modify the software.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Also, for each author's protection and ours, we want to make certain\n");
    (void) fprintf( fp, "that everyone understands that there is no warranty for this free\n");
    (void) fprintf( fp, "software. If the software is modified by someone else and passed on, we\n");
    (void) fprintf( fp, "want its recipients to know that what they have is not the original, so\n");
    (void) fprintf( fp, "that any problems introduced by others will not reflect on the original\n");
    (void) fprintf( fp, "authors' reputations.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Finally, any free program is threatened constantly by software patents.\n");
    (void) fprintf( fp, "We wish to avoid the danger that redistributors of a free program will\n");
    (void) fprintf( fp, "individually obtain patent licenses, in effect making the program\n");
    (void) fprintf( fp, "proprietary. To prevent this, we have made it clear that any patent must\n");
    (void) fprintf( fp, "be licensed for everyone's free use or not licensed at all.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "The precise terms and conditions for copying, distribution and\n");
    (void) fprintf( fp, "modification follow.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "0. This License applies to any program or other work which contains a\n");
    (void) fprintf( fp, "notice placed by the copyright holder saying it may be distributed under\n");
    (void) fprintf( fp, "the terms of this General Public License. The \"Program\", below, refers\n");
    (void) fprintf( fp, "to any such program or work, and a \"work based on the Program\" means\n");
    (void) fprintf( fp, "either the Program or any derivative work under copyright law: that is\n");
    (void) fprintf( fp, "to say, a work containing the Program or a portion of it, either\n");
    (void) fprintf( fp, "verbatim or with modifications and/or translated into another language.\n");
    (void) fprintf( fp, "(Hereinafter, translation is included without limitation in the term\n");
    (void) fprintf( fp, "\"modification\".) Each licensee is addressed as \"you\".\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Activities other than copying, distribution and modification are not\n");
    (void) fprintf( fp, "covered by this License; they are outside its scope. The act of running\n");
    (void) fprintf( fp, "the Program is not restricted, and the output from the Program is\n");
    (void) fprintf( fp, "covered only if its contents constitute a work based on the Program\n");
    (void) fprintf( fp, "(independent of having been made by running the Program). Whether that\n");
    (void) fprintf( fp, "is true depends on what the Program does.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "1. You may copy and distribute verbatim copies of the Program's source\n");
    (void) fprintf( fp, "code as you receive it, in any medium, provided that you conspicuously\n");
    (void) fprintf( fp, "and appropriately publish on each copy an appropriate copyright notice\n");
    (void) fprintf( fp, "and disclaimer of warranty; keep intact all the notices that refer to\n");
    (void) fprintf( fp, "this License and to the absence of any warranty; and give any other\n");
    (void) fprintf( fp, "recipients of the Program a copy of this License along with the Program.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "You may charge a fee for the physical act of transferring a copy, and\n");
    (void) fprintf( fp, "you may at your option offer warranty protection in exchange for a fee.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "2. You may modify your copy or copies of the Program or any portion of\n");
    (void) fprintf( fp, "it, thus forming a work based on the Program, and copy and distribute\n");
    (void) fprintf( fp, "such modifications or work under the terms of Section 1 above, provided\n");
    (void) fprintf( fp, "that you also meet all of these conditions:\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "a) You must cause the modified files to carry prominent notices stating\n");
    (void) fprintf( fp, "that you changed the files and the date of any change.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "b) You must cause any work that you distribute or publish, that in whole\n");
    (void) fprintf( fp, "or in part contains or is derived from the Program or any part thereof,\n");
    (void) fprintf( fp, "to be licensed as a whole at no charge to all third parties under the\n");
    (void) fprintf( fp, "terms of this License.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "c) If the modified program normally reads commands interactively\n");
    (void) fprintf( fp, "when run, you must cause it, when started running for such\n");
    (void) fprintf( fp, "interactive use in the most ordinary way, to print or display an\n");
    (void) fprintf( fp, "announcement including an appropriate copyright notice and a notice that\n");
    (void) fprintf( fp, "there is no warranty (or else, saying that you provide a warranty) and\n");
    (void) fprintf( fp, "that users may redistribute the program under these conditions, and\n");
    (void) fprintf( fp, "telling the user how to view a copy of this License. (Exception: if the\n");
    (void) fprintf( fp, "Program itself is interactive but does not normally print such an\n");
    (void) fprintf( fp, "announcement, your work based on the Program is not required to print an\n");
    (void) fprintf( fp, "announcement.)\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "These requirements apply to the modified work as a whole.\n");
    (void) fprintf( fp, "If identifiable sections of that work are not derived from the Program,\n");
    (void) fprintf( fp, "and can be reasonably considered independent and separate works in\n");
    (void) fprintf( fp, "themselves, then this License, and its terms, do not apply to those\n");
    (void) fprintf( fp, "sections when you distribute them as separate works. But when you\n");
    (void) fprintf( fp, "distribute the same sections as part of a whole which is a work based on\n");
    (void) fprintf( fp, "the Program, the distribution of the whole must be on the terms of this\n");
    (void) fprintf( fp, "License, whose permissions for other licensees extend to the entire\n");
    (void) fprintf( fp, "whole, and thus to each and every part regardless of who wrote it.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Thus, it is not the intent of this section to claim rights or contest\n");
    (void) fprintf( fp, "your rights to work written entirely by you; rather, the intent is to\n");
    (void) fprintf( fp, "exercise the right to control the distribution of derivative or\n");
    (void) fprintf( fp, "collective works based on the Program.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "In addition, mere aggregation of another work not based on the Program\n");
    (void) fprintf( fp, "with the Program (or with a work based on the Program) on a volume of a\n");
    (void) fprintf( fp, "storage or distribution medium does not bring the other work under the\n");
    (void) fprintf( fp, "scope of this License.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "3. You may copy and distribute the Program (or a work based on it, under\n");
    (void) fprintf( fp, "Section 2) in object code or executable form under the terms of Sections\n");
    (void) fprintf( fp, "1 and 2 above provided that you also do one of the following:\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "a) Accompany it with the complete corresponding machine-readable source\n");
    (void) fprintf( fp, "code, which must be distributed under the terms of Sections 1 and 2\n");
    (void) fprintf( fp, "above on a medium customarily used for software interchange; or,\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "b) above on a medium customarily used for software interchange; or, b)\n");
    (void) fprintf( fp, "Accompany it with a written offer, valid for at least three years, to\n");
    (void) fprintf( fp, "give any third party, for a charge no more than your cost of physically\n");
    (void) fprintf( fp, "performing source distribution, a complete machine-readable copy of the\n");
    (void) fprintf( fp, "corresponding source code, to be distributed under the terms of Sections\n");
    (void) fprintf( fp, "1 and 2 above on a medium customarily used for software interchange; or,\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "c) Accompany it with the information you received as to the offer to\n");
    (void) fprintf( fp, "distribute corresponding source code. (This alternative is allowed only\n");
    (void) fprintf( fp, "for noncommercial distribution and only if you received the program in\n");
    (void) fprintf( fp, "object code or executable form with such an offer, in accord with\n");
    (void) fprintf( fp, "Subsection b above.)\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "The source code for a work means the preferred form of the work\n");
    (void) fprintf( fp, "for making modifications to it. For an executable work, complete\n");
    (void) fprintf( fp, "source code means all the source code for all modules it contains,\n");
    (void) fprintf( fp, "plus any associated interface definition files, plus the scripts\n");
    (void) fprintf( fp, "used to control compilation and installation of the executable.\n");
    (void) fprintf( fp, "However, as a special exception, the source code distributed need not\n");
    (void) fprintf( fp, "include anything that is normally distributed (in either source or\n");
    (void) fprintf( fp, "binary form) with the major components (compiler, kernel, and so on) of\n");
    (void) fprintf( fp, "the operating system on which the executable runs, unless that component\n");
    (void) fprintf( fp, "itself accompanies the executable.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "If distribution of executable or object code is made by offering access\n");
    (void) fprintf( fp, "to copy from a designated place, then offering equivalent access to copy\n");
    (void) fprintf( fp, "the source code from the same place counts as distribution of the source\n");
    (void) fprintf( fp, "code, even though third parties are not compelled to copy the source\n");
    (void) fprintf( fp, "along with the object code.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "4. You may not copy, modify, sublicense, or distribute the Program\n");
    (void) fprintf( fp, "except as expressly provided under this License. Any attempt otherwise\n");
    (void) fprintf( fp, "to copy, modify, sublicense or distribute the Program is void, and will\n");
    (void) fprintf( fp, "automatically terminate your rights under this License. However, parties\n");
    (void) fprintf( fp, "who have received copies, or rights, from you under this License will\n");
    (void) fprintf( fp, "not have their licenses terminated so long as such parties remain in\n");
    (void) fprintf( fp, "full compliance.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "5. You are not required to accept this License, since you have not\n");
    (void) fprintf( fp, "signed it. However, nothing else grants you permission to modify or\n");
    (void) fprintf( fp, "distribute the Program or its derivative works. These actions are\n");
    (void) fprintf( fp, "prohibited by law if you do not accept this License. Therefore, by\n");
    (void) fprintf( fp, "modifying or distributing the Program (or any work based on the\n");
    (void) fprintf( fp, "Program), you indicate your acceptance of this License to do so, and all\n");
    (void) fprintf( fp, "its terms and conditions for copying, distributing or modifying the\n");
    (void) fprintf( fp, "Program or works based on it.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "6. Each time you redistribute the Program (or any work based on the\n");
    (void) fprintf( fp, "Program), the recipient automatically receives a license from the\n");
    (void) fprintf( fp, "original licensor to copy, distribute or modify the Program subject to\n");
    (void) fprintf( fp, "these terms and conditions. You may not impose any further restrictions\n");
    (void) fprintf( fp, "on the recipients' exercise of the rights granted herein. You are not\n");
    (void) fprintf( fp, "responsible for enforcing compliance by third parties to this License.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "7. If, as a consequence of a court judgment or allegation of patent\n");
    (void) fprintf( fp, "infringement or for any other reason (not limited to patent issues),\n");
    (void) fprintf( fp, "conditions are imposed on you (whether by court order, agreement or\n");
    (void) fprintf( fp, "otherwise) that contradict the conditions of this License, they do not\n");
    (void) fprintf( fp, "excuse you from the conditions of this License. If you cannot distribute\n");
    (void) fprintf( fp, "so as to satisfy simultaneously your obligations under this License and\n");
    (void) fprintf( fp, "any other pertinent obligations, then as a consequence you may not\n");
    (void) fprintf( fp, "distribute the Program at all. For example, if a patent license would\n");
    (void) fprintf( fp, "not permit royalty-free redistribution of the Program by all those who\n");
    (void) fprintf( fp, "receive copies directly or indirectly through you, then the only way you\n");
    (void) fprintf( fp, "could satisfy both it and this License would be to refrain entirely from\n");
    (void) fprintf( fp, "distribution of the Program.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "If any portion of this section is held invalid or unenforceable under\n");
    (void) fprintf( fp, "any particular circumstance, the balance of the section is intended to\n");
    (void) fprintf( fp, "apply and the section as a whole is intended to apply in other\n");
    (void) fprintf( fp, "circumstances.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "It is not the purpose of this section to induce you to infringe any\n");
    (void) fprintf( fp, "patents or other property right claims or to contest validity of any\n");
    (void) fprintf( fp, "such claims; this section has the sole purpose of protecting the\n");
    (void) fprintf( fp, "integrity of the free software distribution system, which is implemented\n");
    (void) fprintf( fp, "by public license practices. Many people have made generous\n");
    (void) fprintf( fp, "contributions to the wide range of software distributed through that\n");
    (void) fprintf( fp, "system in reliance on consistent application of that system; it is up to\n");
    (void) fprintf( fp, "the author/donor to decide if he or she is willing to distribute\n");
    (void) fprintf( fp, "software through any other system and a licensee cannot impose that\n");
    (void) fprintf( fp, "choice.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "This section is intended to make thoroughly clear what is believed to be\n");
    (void) fprintf( fp, "a consequence of the rest of this License.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "8. If the distribution and/or use of the Program is restricted in\n");
    (void) fprintf( fp, "certain countries either by patents or by copyrighted interfaces, the\n");
    (void) fprintf( fp, "original copyright holder who places the Program under this License may\n");
    (void) fprintf( fp, "add an explicit geographical distribution limitation excluding those\n");
    (void) fprintf( fp, "countries, so that distribution is permitted only in or among countries\n");
    (void) fprintf( fp, "not thus excluded. In such case, this License incorporates the\n");
    (void) fprintf( fp, "limitation as if written in the body of this License.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "9. The Free Software Foundation may publish revised and/or new versions\n");
    (void) fprintf( fp, "of the General Public License from time to time. Such new versions will\n");
    (void) fprintf( fp, "be similar in spirit to the present version, but may differ in detail to\n");
    (void) fprintf( fp, "address new problems or concerns.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "Each version is given a distinguishing version number. If the Program\n");
    (void) fprintf( fp, "specifies a version number of this License which applies to it and \"any\n");
    (void) fprintf( fp, "later version\", you have the option of following the terms and\n");
    (void) fprintf( fp, "conditions either of that version or of any later version published by\n");
    (void) fprintf( fp, "the Free Software Foundation. If the Program does not specify a version\n");
    (void) fprintf( fp, "number of this License, you may choose any version ever published by the\n");
    (void) fprintf( fp, "Free Software Foundation.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "10. If you wish to incorporate parts of the Program into other free\n");
    (void) fprintf( fp, "programs whose distribution conditions are different, write to the\n");
    (void) fprintf( fp, "author to ask for permission. For software which is copyrighted by the\n");
    (void) fprintf( fp, "Free Software Foundation, write to the Free Software Foundation; we\n");
    (void) fprintf( fp, "sometimes make exceptions for this. Our decision will be guided by the\n");
    (void) fprintf( fp, "two goals of preserving the free status of all derivatives of our free\n");
    (void) fprintf( fp, "software and of promoting the sharing and reuse of software generally.\n\n\n");
// GNU END   (see maintenance script update_license_de-GNU)

    return;
}

/*----------------------------------------------------------------------------*/

void show_warranty( FILE *const fp )
{

    (void) fprintf( fp, "NO WARRANTY\n");
    (void) fprintf( fp, "\n");
// GNU BEGIN   (see maintenance script update_license_de-GNU)
    (void) fprintf( fp, "11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY\n");
    (void) fprintf( fp, " THERE IS NO WARRANTY\n");
    (void) fprintf( fp, "FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN\n");
    (void) fprintf( fp, "OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES\n");
    (void) fprintf( fp, "PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER\n");
    (void) fprintf( fp, "EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n");
    (void) fprintf( fp, "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE\n");
    (void) fprintf( fp, "ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH\n");
    (void) fprintf( fp, "YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL\n");
    (void) fprintf( fp, "NECESSARY SERVICING, REPAIR OR CORRECTION.\n");
    (void) fprintf( fp, "\n");
    (void) fprintf( fp, "12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN\n");
    (void) fprintf( fp, "WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY\n");
    (void) fprintf( fp, "AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR\n");
    (void) fprintf( fp, "DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL\n");
    (void) fprintf( fp, "DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM\n");
    (void) fprintf( fp, "(INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED\n");
    (void) fprintf( fp, "INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF\n");
    (void) fprintf( fp, "THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR\n");
    (void) fprintf( fp, "OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.\n\n\n");

// GNU END   (see maintenance script update_license_de-GNU)
    return;
}

/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
