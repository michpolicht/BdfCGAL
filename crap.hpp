#ifndef STUPID_SRC_CMDPARSER_HPP
#define STUPID_SRC_CMDPARSER_HPP

#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <memory>

// C++RAP - C++ Recursive Argument Processor
namespace crap {

class Exception:
        public std::runtime_error
{
	public:
		explicit Exception(const std::string & what);
};

class ExcessiveCmdException:
        public Exception
{
	public:
		explicit ExcessiveCmdException(const std::string & what);
};

class ArgAlreadySetException:
        public Exception
{
	public:
		explicit ArgAlreadySetException(const std::string & what);
};

class ArgRequiresValueException:
        public Exception
{
	public:
		explicit ArgRequiresValueException(const std::string & what);
};

class UnrecognizedArgException:
        public Exception
{
	public:
		explicit UnrecognizedArgException(const std::string & what, int argNum);

		int argNum() const;

	private:
		int m_argNum;
};

class MissingArgException:
        public Exception
{
	public:
		explicit MissingArgException(const std::string & what);
};

class Arg
{
	friend class Parser;
	friend class ArgGroup;

	public:
	    bool isSet() const;

		const std::string & help() const;

		void setHelp(const std::string & help);

		bool required() const;

		void setRequired(bool required);

	protected:
		explicit Arg(const std::string & help);

		/**
		 * Mark argument as being set.
		 */
		void markSet(const std::string & argName);

		virtual int match(char ** argv, int argc) = 0;

		virtual std::string synopsis() const = 0;

		virtual std::string options() const = 0;

		virtual std::string description() const = 0;

	private:
		std::string m_help;
		bool m_required;
		bool m_set;
};

class ValueArg:
    public Arg
{
	friend class Parser;
	friend class ArgGroup;

	public:
	    explicit ValueArg(const std::string & valueName, const std::string & help = "");

	    const std::string & value() const;

		const std::string & valueName() const;

		ValueArg & setValueName(const std::string & valueName);

		const std::string & defaultValue() const;

		ValueArg & setDefaultValue(const std::string & val);

	protected:
		int match(char ** argv, int argc) override;

		std::string synopsis() const override;

		std::string options() const override;

		std::string description() const override;

		void setValue(const std::string & value);

	private:
		std::string m_valueName;
		std::string m_value;
		std::string m_defaultValue;
};


class KeyArg:
    public Arg
{
	friend class Parser;
	friend class ArgGroup;

	public:
	    typedef std::vector<std::string> AliasesContainer;

	    explicit KeyArg(const std::string & name, const std::string & help = "");

	    const std::string & name() const;

		const AliasesContainer & aliases() const;

		KeyArg & addAlias(const std::string & alias);

	protected:
		int match(char ** argv, int argc) override;

		std::string synopsis() const override;

		std::string options() const override;

		std::string description() const override;

		char gluableChar() const;

	private:
		char m_gluableChar;
		AliasesContainer m_aliases;
};

class KeyValueArg:
    public Arg
{
	friend class Parser;
	friend class ArgGroup;

	public:
	    typedef std::vector<std::string> AliasesContainer;

	    KeyValueArg(const std::string & name, const std::string & valueName, const std::string & help = "");

		const std::string & name() const;

		const AliasesContainer & aliases() const;

		KeyValueArg & addAlias(const std::string & alias);

		const std::string & value() const;

		const std::string & valueName() const;

		KeyValueArg & setValueName(const std::string & valueName);

		const std::string & defaultValue() const;

		KeyValueArg & setDefaultValue(const std::string & val);

	protected:
		int match(char ** argv, int argc) override;

		std::string synopsis() const override;

		std::string options() const override;

		std::string description() const override;

		void setValue(const std::string & value);

	private:
		AliasesContainer m_aliases;
		std::string m_valueName;
		std::string m_value;
		std::string m_defaultValue;
};

class Parser;

/**
 * Argument group.
 *
 * There are three types of arguments: key-value arguments (KeyValueArg), value-only arguments (ValueArg) and key-only arguments or flags (KeyArg).
 *
 * Additionaly there are two categories of arguments: commands and attributes. Main differences are as follows.
 *		- Only a single non-required command can be present (otherwise ExcessiveCmdException will be thrown).
 *		- Attributes act as terminals, while commands can have their own groups of sub-arguments.
 *		- Attribute types are distinguishable, which enables type-specific features to be used with them (such as KeyArg gluing).
 *		.
 */
class ArgGroup
{
	friend class Parser;

	public:
	    ArgGroup(const std::string & name = "");

		void setName(const std::string & name);

		std::string name() const;

		void setOptionRequired(bool optionRequired);

		bool optionRequired() const;

		ArgGroup & addAttr(ValueArg * arg);

		ArgGroup & addAttr(KeyArg * arg);

		ArgGroup & addAttr(KeyValueArg * arg);

		Parser * addCmd(Arg * cmd);

	protected:
		typedef std::vector<std::unique_ptr<Parser>> ParsersContainer;
		typedef std::vector<ValueArg *> ValueAttrsContainer;
		typedef std::vector<KeyArg *> KeyAttrsContainer;
		typedef std::vector<KeyValueArg *> KeyValueAttrsContainer;

		void markOptionSet(Arg * cmd);

		Arg * optionSet() const;

		ParsersContainer & parsers();

		ValueAttrsContainer & valueAttrs();

		KeyAttrsContainer & keyAttrs();

		KeyValueAttrsContainer & keyValueAttrs();

		bool gluedKeyArgs(const char * arg) const;

		std::string optionalCmdsSynopsis() const;

		std::string synopsis(std::map<const void *, std::string> & synopsisLines) const;

		std::string description(std::map<const void *, std::string> & descriptionParagraphs) const;

	private:
		std::string m_name;
		bool m_optionRequired;
		Arg * m_optionSet;
		ParsersContainer m_parsers;
		ValueAttrsContainer m_valueAttrs;
		KeyAttrsContainer m_keyAttrs;
		KeyValueAttrsContainer m_keyValueAttrs;
};

/**
 * Arguments parser. Each parser is associated with one command argument. To proceed with parsing the command argument must
 * match the argument passed to the program. Root parser should define program name (argv[0]) as its command argument to be
 * able to process the rest of arguments. The command argument can be set by setCmd() function or it can be passed to the constructor.
 *
 * Command argument can have any number of sub-arguments. Sub-arguments are organized into groups, which are defined by ArgGroup class.
 * For convenience parser provides a default group of sub-arguments. Additional argument groups can be added to the parser. During
 * parsing a parser will try to match arguments defined within the groups. If groups contain command arguments, they will be processed
 * recursively.
 */
class Parser
{
	friend class ArgGroup;

	public:
	    static constexpr char GLUE_CHAR = '-';

		Parser(Arg * cmdArg);

		void setOptionRequired(bool cmdRequired);

		bool optionRequired() const;

		Parser & addAttr(ValueArg * arg);

		Parser & addAttr(KeyArg * arg);

		Parser & addAttr(KeyValueArg * arg);

		Parser & addArgGroup(ArgGroup * argGroup);

		ArgGroup * group(std::size_t index);

		Parser * addSubCmd(Arg * cmd);

		void setHeader(const std::string & header);

		void setFooter(const std::string & footer);

		Arg * cmd() const;

		void setCmd(Arg * cmdArg);

		void printSynopsis(std::ostream & stream = std::cout) const;

		void printDescription(std::ostream & stream = std::cout) const;

		void printHelp(std::ostream & stream = std::cout) const;

		int parse(int argc, char * argv[]);

	protected:
		std::string synopsis(std::map<const void *, std::string> & synopsisLines) const;

		std::string description(std::map<const void *, std::string> & descriptionParagraphs) const;

	private:
		typedef std::vector<ArgGroup *> ArgGroupsContainer;

		Arg * m_cmd;
		ArgGroupsContainer m_argGroups;
		ArgGroup m_defaultGroup;
		std::string m_header;
		std::string m_footer;
};

inline
Exception::Exception(const std::string & what):
    std::runtime_error(what)
{
}

inline
ExcessiveCmdException::ExcessiveCmdException(const std::string & what):
    Exception(what)
{
}

inline
ArgAlreadySetException::ArgAlreadySetException(const std::string & what):
    Exception(what)
{
}

inline
ArgRequiresValueException::ArgRequiresValueException(const std::string & what):
    Exception(what)
{
}

inline
UnrecognizedArgException::UnrecognizedArgException(const std::string & what, int argNum):
    Exception(what),
    m_argNum(argNum)
{
}

inline
int UnrecognizedArgException::argNum() const
{
	return m_argNum;
}

inline
MissingArgException::MissingArgException(const std::string & what):
    Exception(what)
{
}

inline
bool Arg::isSet() const
{
	return m_set;
}

inline
const std::string & Arg::help() const
{
	return m_help;
}

inline
void Arg::setHelp(const std::string & help)
{
	m_help = help;
}

inline
bool Arg::required() const
{
	return m_required;
}

inline
void Arg::setRequired(bool required)
{
	m_required = required;
}

inline
Arg::Arg(const std::string & help):
    m_help(help),
    m_required(false),
    m_set(false)
{
}

inline
void Arg::markSet(const std::string & argName)
{
	if (isSet())
		throw ArgAlreadySetException(std::string() + "Command line argument \"" + argName + "\" has been already set.");
	m_set = true;
}



inline
ValueArg::ValueArg(const std::string & valueName, const std::string & help):
    Arg(help),
    m_valueName(valueName),
    m_value(),
    m_defaultValue()
{
}

inline
const std::string & ValueArg::value() const
{
	if (m_value.empty())
		return m_defaultValue;
	return m_value;
}

inline
const std::string & ValueArg::valueName() const
{
	return m_valueName;
}

inline
ValueArg & ValueArg::setValueName(const std::string & valueName)
{
	m_valueName = valueName;
	return *this;
}

inline
const std::string & ValueArg::defaultValue() const
{
	return m_defaultValue;
}

inline
ValueArg & ValueArg::setDefaultValue(const std::string & val)
{
	m_defaultValue = val;
	return *this;
}

inline
int ValueArg::match(char ** argv, int )
{
	if (!isSet()) {
		setValue(argv[0]);
		return 1;
	}
	return 0;
}

inline
std::string ValueArg::synopsis() const
{
	return std::string() + "<" + valueName() + ">";
}

inline
std::string ValueArg::options() const
{
	if (required())
		return std::string() + " <" + valueName() + "> ";
	else
		return std::string() + "[ <" + valueName() + "> ]";
}

inline
std::string ValueArg::description() const
{
	return std::string(help()).append(" Default value: \"").append(defaultValue()).append("\".");
}

inline
void ValueArg::setValue(const std::string & value)
{
	m_value = value;
	markSet(valueName());
}


inline
KeyArg::KeyArg(const std::string & name, const std::string & help):
    Arg(help),
    m_gluableChar('\0')
{
	addAlias(name);
}

inline
const std::string & KeyArg::name() const
{
	return m_aliases[0];
}

inline
const KeyArg::AliasesContainer & KeyArg::aliases() const
{
	return m_aliases;
}

inline
KeyArg & KeyArg::addAlias(const std::string & alias)
{
	// Check if alias can be glued.
	if ((alias.length() == 2) && (alias[0] == Parser::GLUE_CHAR))
		m_gluableChar = alias[1];
	m_aliases.push_back(alias);
	return *this;
}

inline
int KeyArg::match(char ** argv, int )
{
	for (AliasesContainer::const_iterator it = m_aliases.begin(); it != m_aliases.end(); ++it) {
		if (*it == argv[0]) {
			markSet(argv[0]);
			return 1;
		}
	}
	return 0;
}

inline
std::string KeyArg::synopsis() const
{
	return name();
}

inline
std::string KeyArg::options() const
{
	std::string result;
	if (!required())
		result += "[";
	for (const auto & alias : aliases())
		result.append(" ").append(alias);
	if (!required())
		result += " ]";
	return result;
}

inline
std::string KeyArg::description() const
{
	return help();
}

inline
char KeyArg::gluableChar() const
{
	return m_gluableChar;
}

inline
KeyValueArg::KeyValueArg(const std::string & name, const std::string & valueName, const std::string & help):
    Arg(help),
    m_aliases{name},
    m_valueName(valueName),
    m_value(),
    m_defaultValue()
{
}

inline
const std::string & KeyValueArg::name() const
{
	return m_aliases[0];
}

inline
const KeyValueArg::AliasesContainer & KeyValueArg::aliases() const
{
	return m_aliases;
}

inline
KeyValueArg & KeyValueArg::addAlias(const std::string & alias)
{
	m_aliases.push_back(alias);
	return *this;
}

inline
const std::string & KeyValueArg::value() const
{
	if (m_value.empty())
		return m_defaultValue;
	return m_value;
}

inline
const std::string & KeyValueArg::valueName() const
{
	return m_valueName;
}

inline
KeyValueArg & KeyValueArg::setValueName(const std::string & valueName)
{
	m_valueName = valueName;
	return *this;
}

inline
const std::string & KeyValueArg::defaultValue() const
{
	return m_defaultValue;
}

inline
KeyValueArg & KeyValueArg::setDefaultValue(const std::string & val)
{
	m_defaultValue = val;
	return *this;
}

inline
void KeyValueArg::setValue(const std::string & value)
{
	m_value = value;
	markSet(name());
}

inline
int KeyValueArg::match(char ** argv, int argc)
{
	std::string arg = argv[0];
	std::string val;

	// Check if argument is in form arg=val.
	std::size_t assignPos = arg.find_first_of('=');
	if (assignPos != std::string::npos) {
		val = arg.substr(assignPos + 1);
		arg = arg.substr(0, assignPos);
	}

	for (AliasesContainer::const_iterator it = m_aliases.begin(); it != m_aliases.end(); ++it) {
		if (*it == arg) {
			if (assignPos == std::string::npos) {
				if (argc > 1) {
					if (argv[1][0] == Parser::GLUE_CHAR)
						throw Exception(std::string("Loose argument value can not start with \"") + Parser::GLUE_CHAR + "\" (hint: use arg=value syntax).");
					else {
						setValue(argv[1]);
						return 2;
					}
				} else
					throw ArgRequiresValueException(std::string() + "Command line argument \"" + argv[0] + "\" requires a value.");
			} else {
				setValue(val);
				return 1;
			}
		}
	}
	return 0;
}

inline
std::string KeyValueArg::synopsis() const
{
	return name() + "=<" + valueName() + ">";
}

inline
std::string KeyValueArg::options() const
{
	std::string result;
	if (!required())
		result += "[";
	for (const auto & alias : aliases())
		result.append(" ").append(alias);
	result.append(" <").append(valueName()).append(">");
	if (!required())
		result += " ]";
	return result;
}

inline
std::string KeyValueArg::description() const
{
	return std::string(help()).append(" Default value: \"").append(defaultValue()).append("\".");
}

ArgGroup::ArgGroup(const std::string & name):
    m_name(name),
    m_optionRequired(false),
    m_optionSet(nullptr)
{
}

inline
void ArgGroup::setName(const std::string & name)
{
	m_name = name;
}

inline
std::string ArgGroup::name() const
{
	return m_name;
}

inline
void ArgGroup::setOptionRequired(bool required)
{
	m_optionRequired = required;
}

inline
bool ArgGroup::optionRequired() const
{
	return m_optionRequired;
}

inline
ArgGroup & ArgGroup::addAttr(ValueArg * arg)
{
	m_valueAttrs.push_back(arg);
	return *this;
}

inline
ArgGroup & ArgGroup::addAttr(KeyArg * arg)
{
	m_keyAttrs.push_back(arg);
	return *this;
}

inline
ArgGroup & ArgGroup::addAttr(KeyValueArg * arg)
{
	m_keyValueAttrs.push_back(arg);
	return *this;
}

inline
Parser * ArgGroup::addCmd(Arg * cmd)
{
	m_parsers.push_back(std::unique_ptr<Parser>(new Parser(cmd)));
	return m_parsers.back().get();
}

inline
void ArgGroup::markOptionSet(Arg * cmd)
{
	m_optionSet = cmd;
}

inline
Arg * ArgGroup::optionSet() const
{
	return m_optionSet;
}

inline
ArgGroup::ParsersContainer & ArgGroup::parsers()
{
	return m_parsers;
}

inline
ArgGroup::ValueAttrsContainer & ArgGroup::valueAttrs()
{
	return m_valueAttrs;
}

inline
ArgGroup::KeyAttrsContainer & ArgGroup::keyAttrs()
{
	return m_keyAttrs;
}

inline
ArgGroup::KeyValueAttrsContainer & ArgGroup::keyValueAttrs()
{
	return m_keyValueAttrs;
}

inline
bool ArgGroup::gluedKeyArgs(const char * arg) const
{
	// Basic checks.
	if (arg[0] != Parser::GLUE_CHAR)
		return false;
	if (std::strlen(arg) == 1)
		return false;

	// To be glued key-only arguments each character in the glued string must match one of the key-only arguments.

	// Skip CmdParser::GLUE_CHAR (i == 0).
	for (std::size_t i = 1; i < std::strlen(arg); i++) {
		KeyAttrsContainer::const_iterator it = m_keyAttrs.begin();
		while (it != m_keyAttrs.end()) {
			if (arg[i] == (*it)->gluableChar())
				break;
			++it;
		}
		// These are not glued key-only arguments as one of the characters does not match any of the key-only args.
		if (it == m_keyAttrs.end())
			return false;
	}
	return true;
}

inline
std::string ArgGroup::optionalCmdsSynopsis() const
{
	std::string result;
	for (ParsersContainer::const_iterator it = m_parsers.begin(); it != m_parsers.end(); ++it)
		if (!(*it)->cmd()->required()) {
			if (!result.empty())
				result += '|';
			result += (*it)->cmd()->synopsis();
		}
	return result;
}

inline
std::string ArgGroup::synopsis(std::map<const void *, std::string> & synopsisLines) const
{
	std::string subParserRequiredSynopsis;
	std::string subParserOptionalSynopsis;
	std::string keyRequiredSynopsis;
	std::string keyOptionalSynopsis;
	std::string keyRequiredGluedSynopsis;
	std::string keyOptionalGluedSynopsis;
	std::string keyValRequiredSynopsis;
	std::string keyValOptionalSynopsis;
	std::string valRequiredSynopsis;
	std::string valOptionalSynopsis;

	for (ParsersContainer::const_iterator it = m_parsers.begin(); it != m_parsers.end(); ++it)
		if ((*it)->cmd()->required())
			subParserRequiredSynopsis.append(" ").append((*it)->synopsis(synopsisLines));
	std::string parsersSynopsisString;
	for (ParsersContainer::const_iterator it = m_parsers.begin(); it != m_parsers.end(); ++it)
		if (!(*it)->cmd()->required()) {
			if (!parsersSynopsisString.empty())
				parsersSynopsisString += '|';
			parsersSynopsisString.append((*it)->synopsis(synopsisLines));
		}
	if (!parsersSynopsisString.empty()) {
		if (optionRequired())
			subParserOptionalSynopsis.append(" ").append(parsersSynopsisString);
		else
			subParserOptionalSynopsis.append(" [").append(parsersSynopsisString).append("]");
	}

	for (KeyAttrsContainer::const_iterator it = m_keyAttrs.begin(); it != m_keyAttrs.end(); ++it) {
		if ((*it)->gluableChar() != '\0') {
			if ((*it)->required()) {
				if (keyRequiredGluedSynopsis.empty())
					keyRequiredGluedSynopsis += " -";
				keyRequiredGluedSynopsis += (*it)->gluableChar();
			} else {
				if (keyOptionalGluedSynopsis.empty())
					keyOptionalGluedSynopsis += " [-";
				keyOptionalGluedSynopsis += (*it)->gluableChar();
			}
		} else if ((*it)->required())
			keyRequiredSynopsis.append(" ").append((*it)->synopsis());
		else
			keyOptionalSynopsis.append(" ").append((*it)->synopsis());
	}
	if (!keyOptionalGluedSynopsis.empty())
		keyOptionalGluedSynopsis += "]";

	for (KeyValueAttrsContainer::const_iterator it = m_keyValueAttrs.begin(); it != m_keyValueAttrs.end(); ++it) {
		if ((*it)->required())
			keyValRequiredSynopsis.append(" ").append((*it)->synopsis());
		else
			keyValOptionalSynopsis.append(" [").append((*it)->synopsis()).append("]");
	}

	for (ValueAttrsContainer::const_iterator it = m_valueAttrs.begin(); it != m_valueAttrs.end(); ++it) {
		if ((*it)->required())
			valRequiredSynopsis.append(" ").append((*it)->synopsis());
		else
			valOptionalSynopsis.append(" [").append((*it)->synopsis()).append("]");
	}

	std::string result;
	if (!subParserRequiredSynopsis.empty())
		subParserOptionalSynopsis[0] = '|';
	result.append(subParserRequiredSynopsis).append(subParserOptionalSynopsis);
	result.append(keyRequiredGluedSynopsis).append(keyOptionalGluedSynopsis).append(keyRequiredSynopsis).append(keyOptionalSynopsis);
	result.append(keyValRequiredSynopsis).append(keyValOptionalSynopsis);
	result.append(valRequiredSynopsis).append(valOptionalSynopsis);
	return result;
}

inline
std::string ArgGroup::description(std::map<const void *, std::string> & descriptionParagraphs) const
{
	std::size_t maxWide = 0;

	std::string requiredDescription;
	std::string optionalDescription;

	// Calculate widths for description formatting.
	for (ParsersContainer::const_iterator it = m_parsers.begin(); it != m_parsers.end(); ++it)
		maxWide = std::max((*it)->cmd()->options().length(), maxWide);
	for (KeyAttrsContainer::const_iterator it = m_keyAttrs.begin(); it != m_keyAttrs.end(); ++it)
		maxWide = std::max((*it)->options().length(), maxWide);
	for (KeyValueAttrsContainer::const_iterator it = m_keyValueAttrs.begin(); it != m_keyValueAttrs.end(); ++it)
		maxWide = std::max((*it)->options().length(), maxWide);
	for (ValueAttrsContainer::const_iterator it = m_valueAttrs.begin(); it != m_valueAttrs.end(); ++it)
		maxWide = std::max((*it)->options().length(), maxWide);

	for (ParsersContainer::const_iterator it = m_parsers.begin(); it != m_parsers.end(); ++it) {
		std::string options = (*it)->cmd()->options();
		std::string * description;
		if ((*it)->cmd()->required())
			description = & requiredDescription;
		else
			description = & optionalDescription;
		description->append(options).append(maxWide - options.length(), ' ').append(" - ").append((*it)->cmd()->description()).append("\n");
	}

	for (KeyAttrsContainer::const_iterator it = m_keyAttrs.begin(); it != m_keyAttrs.end(); ++it) {
		std::string options = (*it)->options();
		std::string * description;
		if ((*it)->required())
			description = & requiredDescription;
		else
			description = & optionalDescription;
		description->append(options).append(maxWide - options.length(), ' ').append(" - ").append((*it)->description()).append("\n");
	}

	for (KeyValueAttrsContainer::const_iterator it = m_keyValueAttrs.begin(); it != m_keyValueAttrs.end(); ++it) {
		std::string options = (*it)->options();
		std::string * description;
		if ((*it)->required())
			description = & requiredDescription;
		else
			description = & optionalDescription;
		description->append(options).append(maxWide - options.length(), ' ').append(" - ").append((*it)->description()).append("\n");
	}

	for (ValueAttrsContainer::const_iterator it = m_valueAttrs.begin(); it != m_valueAttrs.end(); ++it) {
		std::string options = (*it)->options();
		std::string * description;
		if ((*it)->required())
			description = & requiredDescription;
		else
			description = & optionalDescription;
		description->append(options).append(maxWide - options.length(), ' ').append(" - ").append((*it)->description()).append("\n");
	}

	std::string result = requiredDescription + optionalDescription;
	for (ParsersContainer::const_iterator it = m_parsers.begin(); it != m_parsers.end(); ++it)
		result += (*it)->description(descriptionParagraphs);

	return result;
}


inline
Parser::Parser(Arg * cmdArg):
    m_cmd(cmdArg),
    m_argGroups{& m_defaultGroup}
{
}

inline
Parser & Parser::addAttr(ValueArg * arg)
{
	m_defaultGroup.addAttr(arg);
	return *this;
}

inline
Parser & Parser::addAttr(KeyArg * arg)
{
	m_defaultGroup.addAttr(arg);
	return *this;
}

inline
Parser & Parser::addAttr(KeyValueArg * arg)
{
	m_defaultGroup.addAttr(arg);
	return *this;
}


inline
Parser & Parser::addArgGroup(ArgGroup * argGroup)
{
	m_argGroups.push_back(argGroup);
	return *this;
}

inline
ArgGroup * Parser::group(std::size_t index)
{
	return m_argGroups.at(index);
}

inline
Parser * Parser::addSubCmd(Arg * cmd)
{
	return m_defaultGroup.addCmd(cmd);
}

inline
void Parser::setOptionRequired(bool cmdRequired)
{
	m_defaultGroup.setOptionRequired(cmdRequired);
}

inline
bool Parser::optionRequired() const
{
	return m_defaultGroup.optionRequired();
}

inline
void Parser::setHeader(const std::string & header)
{
	m_header = header;
}

inline
void Parser::setFooter(const std::string & footer)
{
	m_footer = footer;
}

inline
void Parser::setCmd(Arg * cmdArg)
{
	m_cmd = cmdArg;
}

inline
Arg * Parser::cmd() const
{
	return m_cmd;
}

inline
void Parser::printSynopsis(std::ostream & stream) const
{
	std::map<const void *, std::string> synopsisLines;
	stream << "Usage: " << synopsis(synopsisLines) << "\n";
	for (auto line = synopsisLines.begin(); line != synopsisLines.end(); ++line)
		stream << "       " << line->second << "\n";
}

inline
void Parser::printDescription(std::ostream & stream) const
{
	stream << m_cmd->description() << "\n";
	std::map<const void *, std::string> descriptionParagraphs;
	stream << description(descriptionParagraphs);
	for (auto paragraph = descriptionParagraphs.begin(); paragraph != descriptionParagraphs.end(); ++paragraph)
		stream << paragraph->second;
}

inline
void Parser::printHelp(std::ostream & stream) const
{
	stream << m_header;
	printSynopsis(stream);
	printDescription(stream);
	stream << m_footer;
}

inline
int Parser::parse(int argc, char * argv[])
{
	int argNum = m_cmd->match(argv, argc);
	if (!argNum)
		throw UnrecognizedArgException(std::string() + "Unrecognized argument \"" + argv[argNum] + "\".", argNum);

	while (argNum < argc) {
		int argAdvance = 0;

		for (ArgGroupsContainer::iterator grIt = m_argGroups.begin(); grIt != m_argGroups.end(); ++grIt) {
			ArgGroup * group = *grIt;

			// Look up subcommands first.
			for (ArgGroup::ParsersContainer::iterator it = group->parsers().begin(); it != group->parsers().end(); ++it) {
				try {
					argAdvance = ((*it)->parse(argc - argNum, argv + argNum));
				} catch (UnrecognizedArgException & e) {
					argAdvance = e.argNum();
				}
				if (argAdvance) {
					if (!(*it)->cmd()->required()) {
						if (group->optionSet())
							throw ExcessiveCmdException(std::string() + "Can not use both: \"" + group->optionSet()->synopsis() + "\" and \"" + (*it)->cmd()->synopsis() + "\" at the same time.");
						group->markOptionSet((*it)->cmd());
					}
					break;
				}
			}

			// Check key-value arguments.
			if (!argAdvance)
				for (ArgGroup::KeyValueAttrsContainer::iterator it = group->keyValueAttrs().begin(); it != group->keyValueAttrs().end(); ++it) {
					argAdvance = ((*it)->match(argv + argNum, argc - argNum));
					if (argAdvance)
						break;
				}

			// Check key-only arguments.
			if (!argAdvance)
				for (ArgGroup::KeyAttrsContainer::iterator it = group->keyAttrs().begin(); it != group->keyAttrs().end(); ++it) {
					argAdvance = (*it)->match(argv + argNum, argc - argNum);
					if (argAdvance)
						break;
				}
			// Check if these are glued key-only arguments.
			if (!argAdvance && group->gluedKeyArgs(argv[argNum])) {
				for (std::size_t i = 1; i < std::strlen(argv[argNum]); i++) {
					char glueArg[] = "- ";
					char * glueArgv[] = {glueArg};
					for (ArgGroup::KeyAttrsContainer::iterator it = group->keyAttrs().begin(); it != group->keyAttrs().end(); ++it) {
						glueArg[1] = argv[argNum][i];
						(*it)->match(glueArgv, argc - argNum);
					}
				}
				argAdvance++;
			}

			// If argument does not start with CmdParser::GLUE_CHAR, then handle value-only arguments as it may be one of them.
			if ((!argAdvance) && (argv[argNum][0] != Parser::GLUE_CHAR))
				for (ArgGroup::ValueAttrsContainer::iterator it = group->valueAttrs().begin(); it != group->valueAttrs().end(); ++it) {
					argAdvance = (*it)->match(argv + argNum, argc - argNum);
					if (argAdvance)
						break;
				}

			if (argAdvance)
				break;
		}
		if (argAdvance)
			argNum += argAdvance;
		else
			throw UnrecognizedArgException(std::string() + "Unrecognized argument \"" + argv[argNum] + "\".", argNum);
	}

	for (ArgGroupsContainer::iterator grIt = m_argGroups.begin(); grIt != m_argGroups.end(); ++grIt) {
		ArgGroup * group = *grIt;

		if (group->optionRequired() && !group->optionSet())
			throw MissingArgException(std::string("One of the following arguments must be present: \"") + group->optionalCmdsSynopsis() + "\".");

		// Check if all required arguments are set.
		for (ArgGroup::ParsersContainer::iterator it = group->parsers().begin(); it != group->parsers().end(); ++it)
			if ((*it)->cmd()->required() && !(*it)->cmd()->isSet())
				throw MissingArgException(std::string("Missing required argument \"") + (*it)->cmd()->synopsis() + "\".");
		for (ArgGroup::KeyAttrsContainer::iterator it = group->keyAttrs().begin(); it != group->keyAttrs().end(); ++it)
			if ((*it)->required() && !(*it)->isSet())
				throw MissingArgException(std::string("Missing required argument \"") + (*it)->synopsis() + "\".");
		for (ArgGroup::KeyValueAttrsContainer::iterator it = group->keyValueAttrs().begin(); it != group->keyValueAttrs().end(); ++it)
			if ((*it)->required() && !(*it)->isSet())
				throw MissingArgException(std::string("Missing required argument \"") + (*it)->synopsis() + "\".");
		for (ArgGroup::ValueAttrsContainer::iterator it = group->valueAttrs().begin(); it != group->valueAttrs().end(); ++it)
			if ((*it)->required() && !(*it)->isSet())
				throw MissingArgException(std::string("Missing required argument \"") + (*it)->synopsis() + "\".");
	}

	return argNum;
}

inline
std::string Parser::synopsis(std::map<const void *, std::string> & synopsisLines) const
{
	std::string result(m_cmd->synopsis());
	for (ArgGroupsContainer::const_iterator it = m_argGroups.begin(); it != m_argGroups.end(); ++it) {
		if ((*it)->name().empty())
			result += (*it)->synopsis(synopsisLines);
		else {
			result.append(" (").append((*it)->name()).append(")");
			if (synopsisLines.find(*it) == synopsisLines.end())
				synopsisLines[*it] = std::string("(") + (*it)->name() + ") :=" + (*it)->synopsis(synopsisLines);
		}
	}
	return result;
}

inline
std::string Parser::description(std::map<const void *, std::string> & descriptionParagraphs) const
{
	std::string description;
	for (ArgGroupsContainer::const_iterator it = m_argGroups.begin(); it != m_argGroups.end(); ++it)
		if ((*it)->name().empty())
			description += (*it)->description(descriptionParagraphs);
		else {
			if (descriptionParagraphs.find(*it) == descriptionParagraphs.end())
				descriptionParagraphs[*it] = std::string("(") + (*it)->name() + "):\n" + (*it)->description(descriptionParagraphs);
		}

	if (!description.empty())
		return m_cmd->synopsis() + " options:\n" + description;
	return description;
}

}

#endif
