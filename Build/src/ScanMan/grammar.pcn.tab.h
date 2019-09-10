/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     TOKCOMPOSITION = 258,
     TOKEND = 259,
     TOKFLOAT = 260,
     TOKINT = 261,
     TOKPHASE = 262,
     TOKSPECIES = 263,
     TOKCRAP = 264,
     TOKTHERMO = 265,
     RIGHT_ONE = 266,
     RIGHT_TWO = 267,
     RIGHT_THREE = 268,
     RIGHT_FOUR = 269
   };
#endif
/* Tokens.  */
#define TOKCOMPOSITION 258
#define TOKEND 259
#define TOKFLOAT 260
#define TOKINT 261
#define TOKPHASE 262
#define TOKSPECIES 263
#define TOKCRAP 264
#define TOKTHERMO 265
#define RIGHT_ONE 266
#define RIGHT_TWO 267
#define RIGHT_THREE 268
#define RIGHT_FOUR 269




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 29 "ParseChemkinThermo.y"
{
	int number_i;
	double number_f;
	char string[1000];
}
/* Line 1529 of yacc.c.  */
#line 83 "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/ScanMan/grammar.pcn.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE pcnlval;

