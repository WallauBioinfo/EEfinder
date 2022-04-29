def print_paper_info():
    print("\n")
    print("Mensagem de Cabeçalho")
    print("\n")

def print_end_info():
    print('\n')
    print("|"+"-"*48+"FINISHING"+"-"*48+"|")
    print("\n")
    print("Thank You for use EEfinder")
    print("Please Cite: Paper info")
    print("\n")


class PaperInfo():
    def print_message(self):
        print_paper_info()
    
    def print_finish(self):
        print_end_info()