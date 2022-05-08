#!/usr/bin/env bash

# --- colors --- #

#export NEWT_COLORS='
#root=gray,gray
#shadow=gray,gray
#window=,gray
#border=white,gray
#textbox=white,red
#button=black,white
#'

# --- functions --- #

main_menu() {

	whiptail --title "Main menu" --yesno --yes-button "Proceder au calcul" --no-button "Sortir" "Bienvenue dans l'assistant de résolution de \
	 l'équation de Shrodinger dans le cas du puits de potentiel. \n \nVous allez devoir modifier le fichier nommé 'Fonction_de_potentiel.c'
	 \n ou le laisser intact pour calculer le cas d'un potentiel nul" 15 80
	if [ $? -eq 0 ] ; then
		gcc -c Fonction_de_potentiel.c -o src/Fonction_de_potentiel.o
		gcc -o src/main src/main.c -lgsl -lgslcblas -lm src/Fonction_de_potentiel.o #-Wall -Wextra
		rm src/Fonction_de_potentiel.o
		./src/main
		mv MyFile.csv resultats
		mv src/main src/executable
		gnuplot src/out.p &
		sleep 0.1
		mv graph.png resultats
		
		pkill feh
		feh resultats/graph.png &
	else
		exit 
	fi
}

main_menu



