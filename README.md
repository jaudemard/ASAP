# ASAP

The **ASAP** (**Acessible Surface Area Project**) propose a quick and easy tool to calculate a protein's accessible surface area to a solvant. The project is conducted in an academic context, for the Université Paris Cité Bioinformatic Master Advanded Programming class.

## Introduction

SUJET : CALCUL DE LA SURFACE ACCESSIBLE AU SOLVANT D’UNE PROTEINE
Objectif : Objectif : Réalisez un programme permettant de calculer la surface accessible au solvant 
(absolue et relative) à partir des coordonnées d’une protéine issue d’un fichier PDB. 
Etapes :
- Extraction des coordonnées de chaque atome 
- A partir de chaque atome contenu caractérisant la protéine, créez un nuage de points uniformément 
sur la surface d’une sphère centrée sur l’atome. La sphère aura pour valeur de rayon le rayon de Van 
der Waals de l’atome + rayon de la sonde (solvant = rayon d’un atome d’oxygène) (l’algorithme de 
Saff et Kuijlaars pourra être utilisé).
- A partir de chaque point de la sphère recherchez s’il existe un point appartenant à une autre sphère 
à une distance inférieure au rayon de la sonde 
- A partir des points sans contact, calculez la surface accessible totale en terme de points puis la 
convertir en Å2
 de chaque résidu de la protéine. Calculez ensuite la surface accessible relative, et 
enfin le pourcentage d’accessibilité au solvant.
- Comparez et évaluez la méthode par rapport à NACCESS

## Requirements

**WARNING about NACCESS**
Download of NACCESS require decryption with key sent by mail, sent by the developper of NACCESS, Simon HUBBARD. Given the deadline of the project and the availability of NACCESS on the University computers its installation is not managed here, either ask for the decryption key yourself **or** use a machine with the software already installed.

## Basic Usage

## Description

## Contribution

## Citations

